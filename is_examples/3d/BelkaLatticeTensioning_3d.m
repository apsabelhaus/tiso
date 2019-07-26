%% BelkaLatticeTensioning_3d.m
% Copyright Andrew P. Sabelhaus 2018-2019

% This script used the tiso libraries to calculate cable tensions for a
% desired pose of the Belka quadruped robot, consisting of a tensegrity
% spine with a set of legs.

%% set up the workspace
clear all;
close all;
clc;

% add the core libraries
addpath( genpath('../../is_core/3d') );
% and the helper functions
addpath( genpath('../../is_core/helper') );
% same for the plotting
addpath( genpath('../../is_plot/3d') );
addpath( genpath('../../is_plot/helper') );
% the trajectories are now in a subfolder
addpath( genpath('is_state_traj_3d') );

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 1;

% Startup message if appropriate
if debugging >= 1
    disp('Starting Belka lattice tensioning via inverse statics optimization...');
end

% minimum cable force density
% qMin = 0.5; % units of N/m, depending on m and g
% Putting a higher pretension sometimes helps give consistent bounds on cables and
% prevents too many with different forces.
% qMin = 5;
qMin = 25;
% qMin = 40;

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 3D
% spine. Nodes are column vectors [x; y; z], but for ease, written as transposed
% here.

% FRAME IN MATLAB:
% Belka length axis: -X. Belka's head is facing the -X axis, and hips/rear
% are facing the +X axis.
% Belka width axis: Y. Right shoulder is +Y.
% Belka height axis: Z. 

%%%%%%%%%%%%%%%%%% TO-DO: UPDATE FROM CAD
% bar_endpoint = 0.5; % meters. 50 cm.
be = 142 * 10^-3; % meters. Is 142 mm, from CAD on 2019-07-23.

% This is for the "vertical" spine
% a = [   0,              0,              0;
%         bar_endpoint,     0,              -bar_endpoint;
%         -bar_endpoint,    0,              -bar_endpoint;
%         0,              bar_endpoint,     bar_endpoint;
%         0,              -bar_endpoint,    bar_endpoint]';

% This is for the "horizontal" spine
% columns are x, y, z
a = [   0,      0,    0;
        -be,    0,    -be;
        -be,    0,    be;
        be,     be,   0;
        be,    -be,   0]';    

% For the shoulders (without vertebra attached), this is their local frame,
% with the origin at the center of the vertebra-facing face of the
% shoulders.

% Ordering:
% center between shoulders (origin shifted back in x),
% CoM,
% left shoulder joint,
% right shoulder joint,
% left foot,
% right foot

sh_w = 178 * 10^-3;
sh_h = 275 * 10^-3;
sh_d = 32 * 10^-3;

% sh_com_x = 39 * 10^-3;
% sh_com_z = 12 * 10^-3;
% 2019-07-24: with the 3D print that Albert is doing, parts not quite
% aligned yet, the CoM with respect to the shoulder enclosure has
sh_com_x = 14 * 10^-3;
sh_com_z = 8 * 10^-3;

sh = [-sh_d,        0,      0;
      -sh_com_x,    0,      -sh_com_z;
      -sh_d,        -sh_w,  0;
      -sh_d,        sh_w,   0;
      -sh_d,        -sh_w,  -sh_h;
      -sh_d,        sh_w,   -sh_h]';

% Connecting the shoulders/hips to the vertebrae.
% For the shoulders, shift the "sh" frame so that the origin is at the
% vertebra center, then concatenate the "a" frame.
sh_v_conn = 160 * 10^-3;
% frame with shoulders and vertebra
% the points are column vectors now, not rows.
sh_v_shift = [sh_v_conn; 0; 0;];
a_shoulder = [sh - sh_v_shift, a];

% For the hips, need to rotate the "shoulders" to turn them into hips, then
% concatenate the hips to the last vertebra.
% rotation is around the Z axis by 180 degrees.

% A shortcut to the rotation matrices is MATLAB's nice makehgtform
% command. However, need to truncate off the 4th row/column since the
% command is used for computer graphics and we don't care about
% homongeneity.
rot_z = makehgtform('zrotate', pi);
rot_z = rot_z(1:3, 1:3);

% rotate then translate the shoulders to become hips.
a_hip = [a, (rot_z * sh) + sh_v_shift];

if debugging >= 2
    a
    a_shoulder
    a_hip
end

% number of rigid bodies. This INCLUDES the shoulder/hip assemblies.
b = 5;
% number of nodes.
% This is the number of vertebrae (e.g. total bodies minus 2 for
% shoulder/hip) times nodes per vertebra, plus number of nodes each for
% hips and shoulders.
n = size(a, 2) * (b-2) + size(a_shoulder, 2) + size(a_hip, 2);

% Configuration matrix. Split into source, sink, and body matrices for each
% vertebra.
% Columns are nodes, rows are members (rods or cables)

% source matrix
Cso = [0 1 0 0 0;
       0 0 1 0 0;
       0 0 0 1 0;
       0 0 0 0 1;
       0 0 0 1 0;
       0 0 0 1 0;
       0 0 0 0 1;
       0 0 0 0 1];

% sink matrix
Csi = [0 -1 0 0 0;
       0 0 -1 0 0;
       0 0 0 -1 0;
       0 0 0 0 -1;
       0 -1 0 0 0;
       0 0 -1 0 0;
       0 -1 0 0 0;
       0 0 -1 0 0];

% vertebra matrix
Crb = [1 -1 0 0 0;
       1 0 -1 0 0;
       1 0 0 -1 0;
       1 0 0 0 -1];

% Also now for the shoulders and hips.
% Frame for shoulders is 11 nodes,
% 1. center between shoulders (origin shifted back in x),
% 2. CoM,
% 3. left shoulder joint,
% 4. right shoulder joint,
% 5. left foot,
% 6. right foot
% 7. vertebra center,
% 8-11. v2; ... ; v5 vertebra nodes

% connect vertebra center to shoulder geometric center.
% connect shoulder geometric center to 
% 1. the left shoulder joint,
% 2. right shoulder joint,
% 3. CoM.
% connect the shoulder joints to left / right feet.

% For the shoulders themselves:
Csh = [1  0  0  0  0  0 -1  0  0  0  0;
       1 -1  0  0  0  0  0  0  0  0  0;
       1  0 -1  0  0  0  0  0  0  0  0;
       1  0  0 -1  0  0  0  0  0  0  0;
       0  0  1  0 -1  0  0  0  0  0  0;
       0  0  0  1  0 -1  0  0  0  0  0];

% Concatenate the vertebra matrix into the bottom-right corner here.
% number of rows will be number of rows in vertebra matrix, number columns
% of zeros will be number of nodes in the shoulders only
Csh_v = [Csh;
         zeros(size(Crb,1), size(sh,2)), Crb];

% For the hips, we want to number the members for the vertebra first, then
% the hips. *NOTE: left and right are switched since mirorred.
% Frame for the hips is 11 nodes,
% 1. vertebra center,
% 2-5. v2; ... ; v5 vertebra nodes
% 6. center between hips (origin shifted back in x),
% 7. CoM,
% 8. right hip joint,
% 9. left hip joint,
% 10. right foot,
% 11. left foot

% For the hips themselves:
Chip = [1  0  0  0  0 -1  0  0  0  0  0;
        0  0  0  0  0  1 -1  0  0  0  0;
        0  0  0  0  0  1  0 -1  0  0  0;
        0  0  0  0  0  1  0  0 -1  0  0;
        0  0  0  0  0  0  0  1  0 -1  0;
        0  0  0  0  0  0  0  0  1  0 -1];

% Same internal ordering works though since nodes have the same
% shoulder/hip frame just rotated
Chip_v = [Crb, zeros(size(Crb,1), size(sh,2));
          Chip];

% The bodies' connectivity matrix will be the shoulders, vertebra(e), hips
% along a "diagonal".
% Total number of members is the sum of all number of rows, with b-2
% vertebrae.
r = ((b-2) * size(Crb,1)) + size(Csh_v,1) + size(Chip_v,1);

% Insert into place
Cr = zeros(r, n);
% shoulders
Cr(1:size(Csh_v, 1), 1:size(Csh_v,2)) = Csh_v;
% the multiple vertebrae
for i=1:(b-2)
    % The indices. Row start/end, column start/end.
    % Starts at (shoulders + num prev vert + 1)
    row_s = size(Csh_v, 1) + (i-1)*size(Crb, 1) + 1;
    % ends start + length of vertebra, but need to subtract 1 for
    % off-by-one error with MATLAB's counting.
    row_e = row_s + size(Crb, 1) - 1;
    % columns same format just along the column dimension
    col_s = size(Csh_v, 2) + (i-1)*size(Crb, 2) + 1;
    col_e = col_s + size(Crb, 2) - 1;
    Cr(row_s:row_e, col_s:col_e) = Crb;
end
% the hips at the end
Cr( (end - size(Chip_v,1) + 1):end, (n - size(Chip_v,2) + 1):n) = Chip_v; 

% Get the full connectivity matrix.

% Since we're going to build it up manually, do blocks for the shoulder/hip
% cable connections.
% We need a "source" matrix for the shoulders and a "sink" matrix for the
% hips. Both are the vertebra frames just with a larger number of columns.
Cso_sh = [zeros(size(Cso,1), size(sh,2)), Cso];
Csi_hip = [Csi, zeros(size(Csi,1), size(sh,2))];

% Total number of columns will always be n, total number of nodes.
Cs_sh = zeros(size(Cso, 1), n);
% insert at the leftmost indices
Cs_sh(:, 1:size(Cso_sh,2)) = Cso_sh;
% and for the vertebra sink
Cs_sh(:, size(Cso_sh,2)+1 : size(Cso_sh,2) + size(Csi,2)) = Csi;
% same for the hips, which will be the last-numbered cables.
Cs_hip = zeros(size(Csi, 1), n);
% the hips are at the end
Cs_hip(:, (n-size(Csi_hip,2)+1):n) = Csi_hip;
% the vertebra is second-to-end
Cs_hip(:, (n - size(Csi_hip,2) - size(Cso,2) + 1) : (n - size(Csi_hip,2)) ) = Cso;

% Number of cables is rows of...
% Cs_sh,
% (b-3) * vertebra source or sink, since one pair per vertebra-to-vertebra
% Cs_hip
s = size(Cs_sh,1) + size(Cs_hip,1) + (b-3)*size(Cso,1);

% insert into place
Cs = zeros(s, n);
% shoulders
Cs(1:size(Cs_sh,1), 1:size(Cs_sh,2)) = Cs_sh;
% the multiple vertebrae
% A nice use of the get_C function is to pass in [] as Crb, only generating
% the cable connectivity matrix block for the vertebra-to-vertebra
% connections.
% That's 
Cs_vert = get_C(Cso, Csi, [], b-2);
% if nonempty (i.e. b-2 > 1),
if all(size(Cs_vert) > 0)
    % add in, starting at the bottom-right corner of Cs_sh.
    row_s_s = size(Cs_sh,1) + 1;
    row_e_s = size(Cs_sh,1) + size(Cs_vert,1);
    % for the columns, they need to start and end where the actual elements
    % start and end, not the full columns, since Cs_sh etc. are already
    % expanded to have n columns.
    % that's the last node of the shoulder, i.e, 
    col_s_s = size(Csh_v,2) + 1;
    col_e_s = size(Csh_v,2) + size(Cs_vert,2);
    Cs(row_s_s:row_e_s, col_s_s:col_e_s) = Cs_vert;
end

% the hips at the end
Cs( (end - size(Cs_hip,1) + 1):end, (n - size(Cs_hip,2) + 1):n) = Cs_hip; 

% Concatenate for the full C matrix.
C = [Cs; Cr];

% Now, eta will be a vector with the number of nodes for each body.
eta = zeros(b, 1);
% shoulders and hips
eta(1) = size(Csh_v, 2);
eta(end) = size(Chip_v, 2);
% all the vertebrae
eta(2:end-1) = size(Crb, 2);

% eta = [size(Csh_v,  2);
%        size(Crb,    2);
%        size(Chip_v, 2)];

if debugging >= 2
    C
end

% gravitational constant
g = 9.81;

% The masses here are a bit complicated.
% On 2019-07-18, the middle vertebrae are very light, 100g.
mVert = 0.1;
% On 2019-07-24, shoulders and hips are (from Solidworks)
mShHip = 3.519;

% mass per nodes according to vertebra body
% eta_vert = size(Crb, 2);
% mVertNode = mVert/eta_vert;

% Most points won't have masses on them. We'll only assign mass to the CoM
% (or center of vertebra, respectively.)
m = zeros(n,1);
% Shoulder CoM at node 2 in shoulder frame
m(2) = mShHip;
% For the hips, at node 7 out of 11, or equivalently, end - 4.
m(end-4) = mShHip;

% The vertebrae centers are at node 1 in the "a" frame, which begins right
% after a_shoulder.
% There are b-2 vertebrae.
for i=1:(b-2)
    % this vertebra's center index is 
    % (number of shoulder nodes) + vertebra center node translated out
    index = size(a_shoulder,2) + (i-1)*size(a,2) + 1;
    m(index) = mVert;
end

% and in vector form,
% m = ones(n, 1) * mVertNode;

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 5, and keep nodes 6-10
% (which is vertebra two.)
% w = [0; 0; 0; 0; 0; 1; 1; 1; 1; 1];

% For a b-body structure that removes the leftmost set of nodes,
% zeros for the first body's nodes, ones for the remainder.
% w = [zeros(eta, 1); 
%      kron(ones(b,1), ones(eta, 1))];

% Including all nodes:
% this will be done with Belka since we're calculating reaction forces at
% the feet
w = ones(n,1);

% IMPORTANT! If chosing to remove nodes, must change 'b' also, or else the
% iso optimization problem will FAIL.

%% Trajectory of positions
             
% FOR THE STRAIGHT-LINE BELKA:

% translate the vertebrae out along one axis
% Was be * 3/4 for the 2D spine, but just be is good in 3D.
translation = [ be;
                0;
                0];

% rotation is just local rotation for each vertebra
% For the vertical spine rotating to horizontal:
% rotation = [0; pi/2; 0];
% For the vertebra that has a manually rotated a frame
rotation = [0; 0; 0];

% get the traj for all vertebrae
% As of 2019-07-23, this generates "b" bodies AFTER the initial one at
% zero, so need b-1 for all bodies included here.
% xi = trajStraightMultiBody_3d(translation, rotation, b-1);

% FOR THE CURVED-BACK BELKA:
% starts at x0, last point xf should = (be * (b-1)) translations
x0 = 0;
xf = be * (b-1);
% total height of curved back, from neutral position: arbitrary. try...
back_h = -0.02;

xi = trajArcX_3d(x0, xf, back_h, b);

if debugging >= 2
    xi
end

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, y, z)
% calculate from position trajectory, for each body.
% Body 1: shoulders
coordSh = getCoord_3d(a_shoulder, xi(1:6), debugging);
% Body 2:b-1 = vertebrae
% coordV = getCoord_3d(a, xi(7:12), debugging);
% indices into xi. 6 states per body, so 2nd body starts at 7
xi_v_start = 7;
xi_v_end = 6 + (b-2) * 6;
coordV = getCoord_3d(a, xi(xi_v_start:xi_v_end), debugging);
% Body 3: hips. it's at the end of xi
coordH = getCoord_3d(a_hip, xi(end-6+1:end), debugging);
% stack all back together. One column is a node, so
coord = [coordSh, coordV, coordH];
% coord = getCoord_3d(a, xi, debugging);

if debugging >= 2
    coord
end

% Split up the coordinates.
% The tiso core routine expects these to be column vectors, so a
% transpose is needed also
x = coord(1, :)';
y = coord(2, :)';
z = coord(3, :)';

% Plot for reference
rad = 0.005;
labelsOn = 1;
plotTensegrity3d(C, x, y, z, s, rad, labelsOn);

% Calculate the external reaction forces at the feet.
% Vector of assigning pinned nodes
pinned = zeros(n,1);

% Feet are nodes (in their local frames)
% front left = 5
% front right = 6
% rear left = 10
% rear right = 11
pinned(5) = 1;
pinned(6) = 1;
% the rear are more easily just (end) and (end-1)
pinned(end) = 1;
pinned(end-1) = 1;

% Enforce that the front legs are same and back legs are same
% (see getFrictionlessRxnSymmetric_3d for description of this matrix)
% two constraints
symmetric = zeros(2, n);
% front legs
symmetric(1, 5) = 1;
symmetric(1, 6) = -1;
% back legs
symmetric(2, end-1) = 1;
symmetric(2, end) = -1;

% [rx, ry, rz] = getFrictionlessRxn_3d(x, y, z, pinned, m, g, debugging);
[rx, ry, rz] = getFrictionlessRxnSymmetric_3d(x, y, z, pinned, m, g, symmetric, debugging);

% Initialize external forces
px = rx;
py = ry;
pz = rz;

% Add the gravitational reaction forces for each mass.
% a slight abuse of MATLAB's notation: this is vector addition, no indices
% needed, since py and m are \in R^n.
pz = pz - m*g;

% % Add the gravitational reaction forces for each mass.
% for i=1:n
%     pz(i) = pz(i) - m(i)*g;
% end

% The objective function weighting matrix, via the helper(s).
% For 0.5 q_s'*R*q_s
R = getObj_2norm(s);

%% Solve the inverse statics problem
% invkin_core_3d(x, y, z, px, py, pz, C, s, q_min, debugging}

% % The struct that contains all mandatory inputs to the ISO problem
mandInputs.x = x;
mandInputs.y = y;
mandInputs.z = z;
mandInputs.px = px;
mandInputs.py = py;
mandInputs.pz = pz;
mandInputs.C = C;
mandInputs.b = b;
mandInputs.s = s;

% The optional inputs
optionalInputs.w = w;
optionalInputs.R = R;
optionalInputs.qMin = qMin;
optionalInputs.debugging = debugging;
optionalInputs.eta = eta;

% Solve
[fOpt, qOpt, len, ~, ~] = rbISO_3d(mandInputs, optionalInputs);

% Plot
% rad = 0.01;
% labelsOn = 1;
% plotTensegrity3d(C, x, y, z, s, rad, labelsOn);

% Print out optimal forces

if (debugging >= 1)
    disp('Optimal cable tensions (N) are:');
    fOpt
    disp('(Final) lengths of each cable (meters) are:');
    len
end

% Save these to a CSV file for ease of use.
filename_prefix = 'BelkaLatticeTensioning_';
saveResults = 0;

if saveResults
    disp('Saving Belka Lattice Tensioning Results...');
    % Record the current time in a string
    startTimeString = datestr(datetime('now'));
    % Remove the colons and spaces from this string so that Windows doesn't complain when this repository is cloned
    % colons become dashes, spaces become underscores. Regular expressions to the rescue!
    startTimeString = regexprep(startTimeString, ':', '-');
    startTimeString = regexprep(startTimeString, ' ', '_');
    % the file name
    fileName = strcat(filename_prefix, startTimeString, '.csv');
    % header
    hdr = {};
    hdr{end+1} = 'Belka Lattice Tensioning via Inverse Statics Optimization';
    hdr{end+1} = 'Timestamp:';
    hdr{end+1} = startTimeString;
    hdr{end+1} = 'Pose of the robot:';
%     hdr{end+1} = 'ARCHED BACK';
%     hdr{end+1} = 'STRAIGHT BACK';
    hdr{end+1} = 'Minimum force density (pretensioning) in N/m:';
    hdr{end+1} = num2str(qMin);
    hdr{end+1} = '';
    hdr{end+1} = 'Cable No.,Force (N),Length (m)';
    % concatenate the data together
    cableNo = (1:s)';
    resultsData = [cableNo, fOpt, len];
    % open the file
    fid = fopen(fileName, 'w');
    % write all the header strings
    for j = 1:size(hdr, 2)
        fprintf(fid, '%s\n', string(hdr(j)));
    end
    % close so we can use dlmwrite
    fclose(fid);
    % dlmwrite defaults to comma as a delimiter, nicely.
    dlmwrite(fileName, resultsData, '-append');
end





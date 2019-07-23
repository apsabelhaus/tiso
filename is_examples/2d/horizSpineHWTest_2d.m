%% horizSpineHWTest_2d.m
% Copyright Andrew P. Sabelhaus 2019

% This script used the tiso libraries to calculate the optimal
% equilibrium inputs for a single vertebra, 2D tensegrity spine: the larger
% one used for the spring/summer 2019 data collection.

% The 'trajectory' script here outputs a sequence of equilibrium rest
% length solutions, for use with control.

% Modified from original ..._trajectory script to accomodate for the small
% geometric differences with the hardware test setup.

% ALL UNITS IN METERS

%% set up the workspace
clear all;
close all;
clc;

% add the core libraries
addpath( genpath('../../is_core/2d') );
% and the helper functions
addpath( genpath('../../is_core/helper') );
% same for the plotting
addpath( genpath('../../is_plot/2d') );
addpath( genpath('../../is_plot/helper') );
% and the input trajectories
addpath( genpath('is_state_traj_2d') );

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 1;

% If appropriate, output a starting message.
if debugging >= 1
    disp('Starting hardware test horizontal spine 2d rigid body inverse statics example...');
end

% minimum cable force density
% qMin = 0; % units of N/m, depending on m and g
qMin = 0.5;
% qMin = 1.5;
% qMin = 1.0;
% qMin = 0.1;

% minimum rest length. Some small constant. In meters.
% uMin = 0.001; 
% uMin = 0.2; % This one fails!
% uMin = 0.1; % 10cm.
uMin = 0.05; % 5cm.

% Local frames. We're going to use the getCoord2d script for both
% "vertebrae", but with the "fixed" one having a different frame to account
% for the slightly different exit points of the cables from the eye hooks.

% As per Drew's sketched frames, the moving vertebra will have four nodes,
% with the local origin at the geometric center of the vertebra, and all the
% mass concentrated at node 1 (the inertial center.)
% This results in the plot being a bit off with respect to the physical
% vertebra - the pointy part of the "Y" will be longer in the visualization
% here than in the physical prototype - but the physics are the same.

% 1 = inertial center
% 2 = back bottom, 3 = back top, 4 = right.
% in cm:
c_x = 2.91;
width = 14.42;
f_x = 16.06;
h = 15.24;

% CONVERTED to meters:
c_x = c_x * 10^-2;
width  = width * 10^-2;
f_x = f_x * 10^-2;
h =  h * 10^-2;

% w_1 = 16.06;
% w_2 = 14.42;
% h = 15.24;
% w_5 = 2.91;

a_free = [ -c_x,        0;
           -width,         -h;
           -width,          h;
            f_x,        0]';
        
% In addition, we'll need to know where the dowel pin holes are on the
% vertebra. This is used to get the initial, untensioned pose.
% d_1 = dowel closer to CoM, d_2 = dowel closer to node 4 ("right")
% the y-coordinate here is 0.
d_1 = 0.82;
d_2 = 8.44;

% CONVERTED to meters:
d_1 = d_1 * 10^-2;
d_2 = d_2 * 10^-2;
        
% We now only have three nodes for the fixed vertebra. It's a triangle.
% These are in the GLOBAL frame!
% A = bottom left cable anchor,
% B = top left cable anchor,
% C = center cable anchor (both the saddle cables exit here)
% These are SIGNED so we can do nice things with them later, unlike the
% free vertebra which is not signed (some dimensions don't make signed
% sense.

% Corrected frame for the hardware test:
A_x = -6.56;
A_y = 0;
B_x = -6.56;
B_y = 30.7;
C_x = 18.42;
C_y = 15.24;

% CONVERTED to meters:
A_x = A_x * 10^-2;
A_y = A_y * 10^-2;
B_x = B_x * 10^-2;
B_y = B_y * 10^-2;
C_x = C_x * 10^-2;
C_y = C_y * 10^-2;

% 5 = bottom left (A), 6 = top left (B), 7 = center (C)
a_fixed = [  A_x,       A_y;
             B_x,       B_y;
             C_x,       C_y]';
         
% As with the local vertebra frame, locations of the dowel pins are needed
% to get the initial untensioned vertebra pose
% E = dowel pin closer to anchors (left),
% G = dowel pin on the right

% Corrected frame for the hardware test:
E_x = 33.66;
E_y = 15.24;
G_x = 41.28;
G_y = 15.24;

% CONVERTED to meters:
E_x = E_x * 10^-2;
E_y = E_y * 10^-2;
G_x = G_x * 10^-2;
G_y = G_y * 10^-2;

if debugging >= 2
    a_free
    a_fixed
end

% number of rigid bodies
%b = 2;
% When removing the anchor nodes, it's like removing one of the bodies:
b = 1;

% Configuration matrix for WHOLE STRUCTURE.

% Full connectivity matrix
% Rows 1-4 are cables
% Rows 5-9 are bars
% Columns 1-3 are anchor nodes
% Columns 4-7 are free vertebra nodes

%    A  B  C  c  bb bt f
%    1  2  3  4  5  6  7  
C = [1  0  0  0  -1 0  0;  %  1, cable 1, horizontal bottom
     0  1  0  0  0  -1 0;  %  2, cable 2, horizontal top
     0  0  1  0  -1 0  0;  %  3, cable 3, saddle bottom
     0  0  1  0  0  -1 0;  %  4, cable 4, saddle top
     1  0  -1 0  0  0  0;  %  5, bar 1, bottom anchor triangle leg
     0  1  -1 0  0  0  0;  %  6, ...
     0  0  0  1 -1  0  0;  %  7, bar 3, free vertebra bar 1
     0  0  0  1  0 -1  0;  %  8, ...
     0  0  0  1  0  0 -1]; %  9, bar 6

% Need to specify number of cables, to split up C.
s = 4;
% r follows directly, it's the remainder number of rows.
r = size(C,1) - s;
% ...because C is \in R^{10 x 8}.

% number of nodes
n = size(C, 2);

if debugging >= 2
    C
end

% gravitational constant
% note that we do -mg later for py
g = 9.81;

% Note here that I've used m_i as per-body not per-node.
% Probably accidentally changed notation w.r.t. T-CST 2018 paper.

% Mass of a single vertebra
% Mass as measured with a scale by Jacob on 2019-05-10 was 655g.
mBody = 0.655;

% 2018-12-6: let's see if we just assign all the mass to the center of mass
% node. 2019-50-13: that's node 4, little "c" Center of Mass.
m = zeros(n,1);
m(4) = mBody;
% since we're throwing out the anchor points in the force balance, leave
% those masses as zero.

% Spring constant for the cables. This isn't used in the optimization for
% force (densities) itself, but for conversion into inputs when saving the
% data.
% Here's how to specify 'the same spring constant for all cables'. In N/m.
% For the hardware test, the Jones Spring Co. #241 has a constant of 1.54
% lb/in, which is 270 N/m.
% One set of new McMaster springs is 4.79 lb/in,
% conversion is 
lbin_in_nm = 175.126835;
kappa_i = 4.79 * lbin_in_nm;

% vectorized:
kappa = ones(s,1) * kappa_i;

% On 2018-12-6, we changed cable 3 and 4 to have a higher spring constant, so it
% has less extension, since we were running into hardware limitations.
%kappa(3) = 8.61 * lbin_in_nm;

% with the 25 lb/in spring:
kappa(3) = 25.4 * lbin_in_nm;
kappa(4) = 25.4 * lbin_in_nm;

%kappa(4) = 8.61 * lbin_in_nm;

%% Spring initial length calculations - for collision checking

% Some dimensions of the springs.
% We need to check for collisions / add this to the formulation as a
% constraint.

% The following are all in meters.

% The initial lengths of each of the springs used, from the manufactuter,
% listed here according to their spring constants.
% NOTE that this is TIP-TO-TIP: put the whole spring in the calipers and
% you get this number. This is *NOT* center to center.
l_init_8 = 1 * 2.54 * 0.01;
l_init_4 = 1 * 2.54 * 0.01;
l_init_25 = 0.75 * 2.54 * 0.01;

% The spring extender length (which is added to spring initial length) is 
% 1 inch CENTER TO CENTER of the loops, but the spring hooks to the end of 
% the outer loop, and it's 0.13 inch from center to inside of the
% extender's hook, so total from node (center on one end) to spring tip is
% 1.13 inch. HOWEVER, since the spring is curved, the tip doesn't sit
% exactly at the edge. That extra delta is 1.1 mm (0.0011 m), subtracted
% away.
extender_length = 1.13 * 2.54 * 0.01 - (0.0011);

% We also need to account for the screw tensioner mechanism.
% There are two in use here:
% one with a single screw point, and 
% one with two screws (one spring, one cable.)

% For both of these, we'll make the following assumption:
% Screw is held against edge of spring hook,
% cable is tightened against screw from the other side.
% So, for the single-screw tensioner, we're just looking at the diameter of
% the screw which needs to be subtracted from spring. We have M3 screws.
% If need more precision later: include thickness of spring wire.
tensioner_single_length = 0.003;

% For the single-screw tensioner, the spring tip is its center, 
% adding an extra few mm from its center to edge. Measured roughly (cm -> m)
% TO-DO: adjust this more if we get noisy data. Account for screw diameter,
% use that instead, since the cable is really looped around the screw, then
% also in addition account for the offset from center (edge??) of spring
% loop.
% adjuster_length = 0.36 * 0.01;
% tensioner_single_length = 0.36 * 0.01;

% For the double-screw tensioner, 
% center-to-center of the two screws SHOULD be 9mm, but Drew measured 1 cm.
% minus the radius of the screw on each side equals a full diameter (3mm), 
% so inside-edge of screw to inside-edge is
tensioner_double_length = 0.007;

%initial_lengths = [l_init_4; l_init_4; l_init_8; l_init_4];

% The "initial lengths" here are then total distance from node center to
% connection point of cable.
% Cables 1 and 2 (top and bottom) use the 4 lb/in springs, and cables 3 and
% 4 (saddles) use the 25 lb/in spring.
% Also, 1,2 use double extenders, 3,4 use single extenders.
initial_lengths = [l_init_4 + extender_length + tensioner_double_length; 
                   l_init_4 + extender_length + tensioner_double_length; 
                   l_init_25 + extender_length + tensioner_single_length;
                   l_init_25 + extender_length + tensioner_single_length];

% So, the total length to subtract from the rest length is
% init_len_offset = initial_lengths + extender_length;
% ...now included in "initial_lengths."
init_len_offset = initial_lengths;

%% Anchors

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 3, and keep nodes 4-7
% (which is the free vertebra.)
w = [0; 0; 0; 1; 1; 1; 1];
% Including all nodes:
%w = ones(n,1);

% IMPORTANT! If chosing to remove nodes, must change 'b' also, or else inv
% stat will FAIL.

%% Trajectory of positions

% all the positions of each rigid body (expressed as their COM positions
% and rotation). that's 3 states: [x; y; \gamma] with
% the angle being an intrinsic rotation.

% We need the initial pose of the robot, PI_0, at which the cables are
% calibrated. This is for the hardware test: at pose PI_0, the cables are
% assumed to be at "perfectly 0 force", e.g. just barely not slack. For
% example, if we fix the spine at some position, and adjust each cable so
% it's just barely no longer slack, then we have a calibration.

% For calculations with the hardware test setup, we need to know
% where the vertebra is pinned into place for calibration. This is
% currently (as of 2019-05-13) where the dowel pins align ("d1" in free
% frame and "E" in fixed frame.) This is because the zero of the free frame
% is at the geometric center, now.

calibratedPosFreeX = E_x - d_1;
calibratedPosFreeY = E_y;
calibratedPosFree = [calibratedPosFreeX;
                            calibratedPosFreeY;
                            0];

% The "fixed" vertebra, the leftmost one, has no translation or rotation.
% To get the full system state, we must specify both.
% 2019-05-13: easier to just say the global frame for fixed.
calibratedCoordAnchors = a_fixed;

% and for the free vertebra
calibratedCoordFree = getCoord_2d(a_free, ...
    calibratedPosFree, debugging);

% We can then use the same trick from the is_core functions to calculate
% the lengths of the cables, via the C matrix. Concatenate the node
% positions for the whole structure:
calibratedCoordX = zeros(n,1);
calibratedCoordY = zeros(n,1);
% The outputs from get_node_coordinates are row vectors:
calibratedCoordX(1:3) = calibratedCoordAnchors(1,:)';
calibratedCoordX(4:7) = calibratedCoordFree(1,:)';
calibratedCoordY(1:3) = calibratedCoordAnchors(2,:)';
calibratedCoordY(4:7) = calibratedCoordFree(2,:)';
% The lengths of each cable member in the two directions are (via the 
% "cable" rows of C),
dx0 = C(1:4,:) * calibratedCoordX;
dy0 = C(1:4,:) * calibratedCoordY;
% so the lengths of each cable are the euclidean norm of each 2-vector,
% re-organize:
D0 = [dx0, dy0];
% the row norm here is then the length.
len0 = vecnorm(D0, 2, 2);

% Now, for the trajectory as the robot moves:
% The trajectory generation function gives back a sequence of states for
% which we require cable inputs.
% Let's do a trajectory that sweeps from 0 to pi/8.
maxSweep = pi/8;
% maxSweep = pi/12;
% maxSweep = pi/16;
% To swing the vertebra down slowly, do a sweep to a negative number.
% maxSweep = -pi/16;
% Now, also include a minimum. This was 0 before. Now, we can start "down"
% somewhere and bend upward.
% minSweep = -pi/16;
% minSweep = -pi/12;
minSweep = -pi/8;
% minSweep = 0;

% The frames are at an initial configuration of zero rotation.
rotation0 = 0;

% and translated outward.

% 2019-05-13: With the new coordinate frame, let's move the vertebra back
% by half the length from its origin to the back nodes, keeping its
% y-position the same. That width dimension is "width", so 
translation0 = [calibratedPosFreeX - width/3;
                 calibratedPosFreeY];


% with a large number of points.
% numPts = 400;
% For doing the hardware test: the motor controller doesn't have great
% resolution. So, do a smaller number of poins.
% numPts = 5;
% numPts = 2;
% numPts = 10;
% numPts = 20;

numPts = 50;
% Used for the data collected June 2019
% numPts = 80;

% For the frames we've chosen, it doesn't make sense to rotate the vertebra
% around the origin: we do not want to sweep out from the corner of the
% test setup (the global frame zero.)
% A good location of the axis of rotation might be halfway between the A,B
% back nodes and the C node. Or...?
rotAxisPt = [ (C_x + A_x); C_y];

% get the trajectory:
[xiAll] = trajBend2b_2d(translation0, rotation0, minSweep, ...
    maxSweep, rotAxisPt, numPts);

%% Calculations for the inputs to the core tiso library

% The nodal coordinates (x, z)
% calculate from position trajectory. 
% We can do these out for each point, along the trajectories.
% initialize the results. There are n nodes (n position vectors.)
x = zeros(n, numPts);
y = zeros(n, numPts);

for i=1:numPts
    % We split up according to rigid body. The fixed body is xiAll(1:3,i)
    coordinatesFixed_i = getCoord_2d(a_fixed, ...
        xiAll(1:3, i), debugging);
    % and for the free vertebra
    coordinatesFree_i = getCoord_2d(a_free, ...
        xiAll(4:6, i), debugging);
    % ...and split them into coordinate-wise vectors per node, per
    % vertebra.
    x(1:3, i) = coordinatesFixed_i(1,:)';
    x(4:7, i) = coordinatesFree_i(1,:)';
    y(1:3, i) = coordinatesFixed_i(2,:)';
    y(4:7, i) = coordinatesFree_i(2,:)';
end

if debugging >= 2
    %coordinates
    x
    y
end

% Reaction forces can be calculated by this helper, which assumes that only
% gravity is present as an external force.
% Initialize results
px = zeros(n, numPts);
py = zeros(n, numPts);

% Add the gravitational reaction forces for each mass.
% a slight abuse of MATLAB's notation: this is vector addition, no indices
% needed, since py and m are \in R^n.
py = py - m*g;


%% Solve the inverse statics problem

% the objective function weight is
% R = getObj_PE(s, lengths, kappa);
R = getObj_2norm(s);

% Solve, over each point.
% Let's use a cell array for Ab and pb, since I don't feel like thinking
% over their sizes right now.
fOpt = zeros(s, numPts);
qOpt = zeros(s, numPts);
% we also need the lengths for calculating u from q*.
len = zeros(s, numPts);
Ab = {};
pb = {};

% finally, the big function call:
for i=1:numPts
    if debugging >= 1
        disp('Iteration:');
        disp(num2str(i));
    end
    % combine parameters into the mandatory and optional arguments to core
    % The struct that contains all mandatory inputs to the ISO problem
    mandInputs_i.x = x(:,i);
    mandInputs_i.y = y(:,i);
    mandInputs_i.px = px(:,i);
    mandInputs_i.py = py(:,i);
    mandInputs_i.C = C;
    mandInputs_i.b = b;
    mandInputs_i.s = s;

    % The optional inputs
    optionalInputs_i.w = w;
    optionalInputs_i.R = R;
    optionalInputs_i.qMin = qMin;
    optionalInputs_i.debugging = debugging;
    % including the uMin constraint
    optionalInputs_i.uMin = uMin;
    optionalInputs_i.kappa = kappa;
    
    % quadprog is inside this function.
    [fOpt(:,i), qOpt(:,i), len(:,i), Ab_i, pb_i] = rbISO_2d(mandInputs_i, ...
        optionalInputs_i);
    
    % and insert this Ab and pb.
    Ab{end+1} = Ab_i;
    pb{end+1} = pb_i;
end

% A quick plot of the cable tensions.
figure; 
hold on;
subplot(4,1,1)
hold on;
title('Cable tensions');
plot(fOpt(1,:))
ylabel('1 (N)');
subplot(4,1,2)
plot(fOpt(2,:))
ylabel('2 (N)');
subplot(4,1,3)
plot(fOpt(3,:))
ylabel('3 (N)');
subplot(4,1,4);
plot(fOpt(4,:));
ylabel('4 (N)');

%% Convert the optimal forces into optimal rest lengths.
% u_i = l_i - (F_i/kappa_i)

% save in a vector
uOpt = zeros(s, numPts);
availableLen = zeros(s, numPts);

% it's more intuitive to iterate for now. At least, we can iterate over
% cables and not over timesteps.
for k=1:s
    % For cable k, divide the row in f_opt by kappa(k)
    uOpt(k, :) = len(k,:) - (fOpt(k,:) ./ kappa(k));
    % But, now include the length offset term. Accounts for the initial
    % spring length, as well as the little extender we had to use.    

    % Instead, collisions occur when u_opt - init_len_off < u_min.
    availableLen(k,:) = uOpt(k,:) - init_len_offset(k);
end

% Check to confirm that these satisfy the desired input constraint.
% A better way to do this would be to set u_min(k) to init_len_offset(k) +
% some small number, since it's not really not u > 0, there's the not-cable
% part of the length that can't be spooled past the anchors.

% Does the cable get smaller than u_min?
disp('Cable collisions: 0 = none, 1 = length shorter than min + offset');
availableLen < uMin


% For use with the hardware example, it's easier to instead define a
% control input that's the amount of "stretch" a cable experiences.
% Note that we also save a difference from Pi_0, the calibrated position.
% Some algebra shows that the control input here is stretch added to the
% delta of absolute cable lengths from Pi_0 to lengths_i.
stretchOpt = zeros(s, numPts);
lengthsDiff = zeros(s, numPts);
stretchOptAdj = zeros(s, numPts);

for k=1:s
    % For cable k, divide the row in f_opt by kappa(k)
    stretchOpt(k, :) = (fOpt(k,:) ./ kappa(k));
    % the length difference is a subtraction from the initial pose
    % note that lengths_0 is transposed as calculated above.
    % We're iterating over cables, so need to index into the initial
    % lengths. Abusing some MATLAb broadcasting here (properly, ones(1,k)
    % multiplied by lengths_0(k).)
    lengthsDiff(k,:) = len0(k)' - len(k,:);
    % and then add to the adjusted stretch,
    % ALSO SCALE TO GET CM, since the microcontroller uses that best, not
    % meters.
    stretchOptAdj(k,:) = (stretchOpt(k,:) + lengthsDiff(k,:)) * 100;
    
    % 2019-05-13: NOTE HERE that positive "lengths_diff" means positive
    % rettraction on the cable length: e.g., lengths_diff(1,1) = 0.107 means
    % that cable 1 is 10.7cm *SHORTER* at time 0 versus calibration.
    % Meaning that the motor should be retracting 10.7cm in addition to
    % the required stretch amount (also positive for retraction.)
end

% A quick plot of the changes in cable lengths.
figure; 
hold on;
subplot(4,1,1)
hold on;
title('Cable Inputs (Amount of Retraction From Calibration)');
plot(stretchOptAdj(1,:))
ylabel('1 (cm)');
subplot(4,1,2)
plot(stretchOptAdj(2,:))
ylabel('2 (cm)');
subplot(4,1,3)
plot(stretchOptAdj(3,:))
ylabel('3 (cm)');
subplot(4,1,4);
plot(stretchOptAdj(4,:));
ylabel('4 (cm)');

%% Plot the structure, for reference.

% This should make it easier to visualize the results.
% Need to specify "how big" we want the bars to be. A good number is
% radius = 0.005; % meters.
% 2019-05-13: The actual vertebra is 4 cm in width, so 
radius = 0.02;
% Need to modify the plotting script to make this more sensible.

% Plotting: do the first position and the last position.
% Basic command:
% plotTensegrity2d(C, x, y, s, radius);
% We want to get the first and last coordinates.
% Reference configuration state:
plotTensegrity2d(C, calibratedCoordX, calibratedCoordY, s, radius);
% Initial position:
plotTensegrity2d(C, x(:,1), y(:,1), s, radius);
% Final position:
plotTensegrity2d(C, x(:,end), y(:,end), s, radius);

% In order to use this date

%% Save the data.

% path to store: ***CHANGE THIS PER-USER***
% for now, use the user's home directory.
% savefile_path = '~/';
% Or for more compatibility with different operating systems,
savefilePath = strcat(pwd, filesep);

% write the actual data
% we used the rigid body reformulation method here, 
nOrB = 1;
%saveISOresult2d(uOpt, xi_all, n, r, n_or_b, savefile_path);

% For the hardware test, we want to use "stretch" not rest length.
% saveISOresult2d(stretchOptAdj, xiAll, n, r, nOrB, savefilePath);
















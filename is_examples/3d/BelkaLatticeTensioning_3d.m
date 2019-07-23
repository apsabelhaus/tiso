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
debugging = 2;

% Startup message if appropriate
if debugging >= 1
    disp('Starting Belka lattice tensioning via inverse statics optimization...');
end

% minimum cable force density
% qMin = 0.5; % units of N/m, depending on m and g
% Putting a higher pretension sometimes helps give consistent bounds on cables and
% prevents too many with different forces.
qMin = 3;

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 3D
% spine. Nodes are column vectors [x; y; z], but for ease, written as transposed
% here.

%%%%%%%%%%%%%%%%%% TO-DO: UPDATE FROM CAD
bar_endpoint = 0.5; % meters. 50 cm.

% This is for the "vertical" spine
% a = [   0,              0,              0;
%         bar_endpoint,     0,              -bar_endpoint;
%         -bar_endpoint,    0,              -bar_endpoint;
%         0,              bar_endpoint,     bar_endpoint;
%         0,              -bar_endpoint,    bar_endpoint]';

% This is for the "horizontal" spine
a = [   0,                0,              0;
        -bar_endpoint,    0,              -bar_endpoint;
        -bar_endpoint,    0,              bar_endpoint;
        bar_endpoint,     bar_endpoint,   0;
        bar_endpoint,    -bar_endpoint,   0]';    
    
if debugging >= 2
    a
end

% number of rigid bodies. Note this is number of MOVING bodies now!
% so we have five vertebrae, 4 moving.
% b= 4;
b = 3;
% number of nodes
n = size(a, 2) * b;

% Configuration matrix. Split into source, sink, and body matrices for each
% vertebra,

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
   
% Get the full connectivity matrix.
% NOTE HERE that though we'll be removing the leftmost vertebra later, we
% still need its connectivity matrix blocks, so number of bodies is
% actually b+1.

C = get_C(Cso, Csi, Crb, b);


% Need to specify number of cables, to split up C.
% Cables per pair of bodies is number in source (or sink, equivalently)
sigma = size(Cso, 1);
% and there is one pair of cables between each moving body, meaning (for
% example) two sets of cables for three bodies, so b-1 here
s = sigma * (b-1);

% r follows directly, it's the remainder number of rows.
r = size(C,1) - s;

% later, we will need number of nodes per body. Could get num col of any of
% the submatrices here
eta = size(Cso, 2);

% number of nodes
n = size(C, 2);

if debugging >= 2
    C
end

% gravitational constant
g = 9.81;

% Mass of a single vertebra
% On 2019-07-18, the middle vertebrae are very light, 100g.
% To-do: estimates from Solidworks about mass of shoulders and hips.
mBody = 0.1;

% mass per nodes according to vertebra body
mNode = mBody/eta;

% and in vector form,
m = ones(n, 1) * mNode;

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
             
% translate the vertebrae out along one axis
translation = [ bar_endpoint * (3/4);
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
xi = trajStraightMultiBody_3d(translation, rotation, b-1);

if debugging >= 2
    xi
end

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, y, z)
% calculate from position trajectory
coord = getCoord_3d(a, xi, debugging);

if debugging >= 2
    coord
end

% Split up the coordinates.
% The tiso core routine expects these to be column vectors, so a
% transpose is needed also
x = coord(1, :)';
y = coord(2, :)';
z = coord(3, :)';

% A test: calculate the external reaction forces if certain nodes are
% pinned. That would be 4 and 5 for the left vertebra, and maybe for
% example these same on the 3rd vertebra, which are 14 and 15
pinned = zeros(n,1);
pinned(4) = 1;
pinned(5) = 1;
pinned(14) = 1;
pinned(15) = 1;

% Enforce that the front legs are same and back legs are same
% (see getFrictionlessRxnSymmetric_3d for description of this matrix)
% two constraints
symmetric = zeros(2, n);
% front legs
symmetric(1, 4) = 1;
symmetric(1, 5) = -1;
% back legs
symmetric(2, 14) = 1;
symmetric(2, 15) = -1;

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

% Solve
[fOpt, qOpt, len, ~, ~] = rbISO_3d(mandInputs, optionalInputs);

% Plot
rad = 0.02;
labelsOn = 1;
plotTensegrity3d(C, x, y, z, s, rad, labelsOn);





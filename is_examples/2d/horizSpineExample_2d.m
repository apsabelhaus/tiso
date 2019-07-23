%% horizSpineIS_example_2d.m
% Copyright Andrew P. Sabelhaus 2019

% This script used the tiso libraries to solve the inverse statics problem
% for a two body (b=2) tensegrity spine, defined in d=2 dimensions. 
% As the term is used here, 'spine'
% refers to a structure with repeated rigid bodies of the same geometry and

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

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 1;

% If appropriate, output a starting message.
if debugging >= 1
    disp('Starting horizontal spine 2d rigid body inverse statics example...');
end

% minimum cable force density
%qMin = 0; % units of N/m, depending on m and g
qMin = 0.5;

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 2D
% spine. Nodes are column vectors [x; z], but for ease, written as transposed
% here.

%barEndpoint = 0.5; % meters. 50 cm.

% For the 2D spine control test, fall 2018, the dimensions of the vertebra
% are as follows:
% barEndpoint = 4 inches
% that's 4 * 2.54 * 0.01 = 0.1016 meters
barEndpoint = 0.1016;

% local frame of nodes
a = [   0,              0;
        barEndpoint,   -barEndpoint;
        -barEndpoint,  -barEndpoint;
        0,              barEndpoint]';
    

% (e.g., vertebra is in an 8x8 inch box.)

% Mass of a single vertebra
% Mass as measured with a scale on 2018-11-18 is about 500g
mBody = 0.495;
    
if debugging >= 2
    a
end

% number of rigid bodies
b = 2;
% When removing the anchor nodes, it's like removing one of the bodies:
% b = 1;

% Configuration matrix for WHOLE STRUCTURE.

% Rows 1-4 are cables
% Rows 5-10 are bars
% Columns 1-4 are bottom vertebra nodes
% Columns 5-8 are top vertebra nodes

%    1  2  3  4  5  6  7  8  
C = [0  1  0  0  0 -1  0  0;  %  1, cable 1
     0  0  1  0  0  0 -1  0;  %  2, ...
     0  0  0  1  0 -1  0  0;  %  3, ...
     0  0  0  1  0  0 -1  0;  %  4, cable 4
     1 -1  0  0  0  0  0  0;  %  5, bar 1
     1  0 -1  0  0  0  0  0;  %  6, ...
     1  0  0 -1  0  0  0  0;  %  7, ...
     0  0  0  0  1 -1  0  0;  %  8, ...
     0  0  0  0  1  0 -1  0;  %  9, ...
     0  0  0  0  1  0  0 -1]; % 10, bar 6

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
g = 9.81;

% mass per nodes according to vertebra body
% 4 nodes per body, so
mNode = mBody/4;

% and in vector form,
m = ones(n, 1) * mNode;

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 4, and keep nodes 5-8
% (which is vertebra two.)
% w = [0; 0; 0; 0; 1; 1; 1; 1];
% Including all nodes:
w = ones(n,1);

% IMPORTANT! If chosing to remove nodes, must change 'b' also, or else the
% iso optimization problem will FAIL.

% We also need to declare which nodes are pinned (with external reaction
% forces) and which are not. 

% ******NOTE: this is NOT supported in the version of tiso as of June 2019.
%       The function for obtaining reaction forces has not been verified.

% We're choosing not to have this be the same as w, since there are some
% nodes to "ignore" for w which are not necessarily built in to the ground
% for "pinned". Example, nodes inside the leftmost vertebra, where we're
% deciding to assume that only the tips of its "Y" are supported.

% So, only nodes 2 and 3 are pinned.
pinned = zeros(n,1);
pinned(2) = 1;
pinned(3) = 1;

%% Trajectory of positions

% all the positions of each rigid body (expressed as their COM positions
% and rotation). that's 3 states: [x; z; \gamma] with
% the angle being an intrinsic rotations.

% use \xi for the system states.
xi = zeros(b * 3, 1);

% for rigid body 1, the fixed one, doesn't move. However, we've defined it
% in its "vertical" state, so it needs to be rotated by 90 degrees around
% the y-axis so the saddle cables can align. Let's rotate it clockwise so
% that node 4 is in +x.
xi(3) = -pi/2;

% for rigid body 2, translate out in the +x direction. Translating by one
% full body length puts the tips exactly in the same plane, so maybe do it
% to 3/4 of that length.
% the length of one vert is 2 * bar_endpoint. 
% x-position is coordinate 1.
% To make things interesting, let's rotate it a small bit, too.
%xi(3) = -pi/2 + pi/16;

% xi(4:6) = [     barEndpoint * (3/2);
%                 0;
%                -pi/2 + pi/16];

% No rotation:
xi(4:6) = [     barEndpoint * (3/2);
                0;
               -pi/2];

if debugging >= 2
    xi
end

%% Calculations for the inputs to the core tiso library

% The nodal coordinates (x, z)
% calculate from position trajectory
coord = getCoord_2d(a, xi, debugging);

if debugging >= 2
    coord
end

% Split up the coordinates.
% The tiso core routine expects these to be column vectors, so a
% transpose is needed also
x = coord(1, :)';
y = coord(2, :)';

% Reaction forces can be calculated by this helper, which assumes that only
% gravity is present as an external force.
[px, py] = getRxn2d(x, y, pinned, m, g, debugging);

% for more details, you can look at commits to the library before Nov. 2018
% where this reaction force/moment balance was written out by hand.

% Since this was just a calculation of the reaction forces, we ALSO need to
% add in the external forces themselves (grav forces) for use as the whole
% inverse statics problem.

% Add the gravitational reaction forces for each mass.
% a slight abuse of MATLAB's notation: this is vector addition, no indices
% needed, since py and m are \in R^n.
py = py + -m*g;

% The objective function weighting matrix, via the helper(s).
% For 0.5 q_s'*R*q_s
% R = getObj_2norm(s);

% For the potential energy objective function,
lengths = getLen_2d(x, y, s, C);
% need to declare some spring constants
lbin_in_nm = 175.126835;
kappa_i = 4.79 * lbin_in_nm;
% vectorized:
kappa = ones(s,1) * kappa_i;

% so the objective function weight is
% R = getObj_PE(s, lengths, kappa);
R = getObj_2norm(s);

%% Solve the inverse statics problem

% The struct that contains all mandatory inputs to the ISO problem
mandInputs.x = x;
mandInputs.y = y;
mandInputs.px = px;
mandInputs.py = py;
mandInputs.C = C;
mandInputs.b = b;
mandInputs.s = s;

% The optional inputs
optionalInputs.w = w;
optionalInputs.R = R;
optionalInputs.qMin = qMin;
optionalInputs.debugging = debugging;

% Solve. Rigid body reformulation, two dimensions.
[fOpt, qOpt, ~, Ab, pb] = rbISO_2d(mandInputs, optionalInputs);


% [fOpt, qOpt, Ab, pb] = rbISO_2d(x, y, px, py, w, C, R, s, b, qMin, debugging);

% Seems correct, intuitively!
% Cable 1 is horizontal, below.
% Cable 2 is horizontal, above.
% Cable 3 is saddle, below.
% Cable 4 is saddle, above.

% It makes sense that cable 2 force > cable 1 force, for an "upward" force
% in the saggital plane, counteracting gravity on the suspended vertebra.
% it then also could make sense that cable 3 has no force: the
% gravitational force on the suspended vertebra gives a clockwise moment,
% which is also the direction of cable 3.
% Cable 4's force is likely counteracting the gravitational moment of the
% suspending vertebra.

% TO-DO: Confirm, by hand/in simulation, that these forces actually keep
% the vertebra in static equilibrium.

% TO-DO: sign check. Are we applying positive tension forces or negative
% tension forces? See if/what solution pops out with the constraint in the
% opposite direction (< c, not > c.)


%% Plot the structure, for reference.

% This should make it easier to visualize the results.
% Need to specify "how big" we want the bars to be. A good number is
%radius = 0.02; % meters.
% For the hardware test in fall 2018, a better number, based on the
% hardware itself, is 1.5 inches, which in meters is
%radius = 0.038;
% but that's too big given how the plotter works right now. Choose
% something arbitrarily smaller.
radius = 0.01;

% Plot.
plotTensegrity2d(C, x, y, s, radius);

















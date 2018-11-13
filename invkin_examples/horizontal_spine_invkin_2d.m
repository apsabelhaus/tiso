%% horizontal_spine_invkin_2d.m
% Copyright Andrew P. Sabelhaus 2018

% This script used the tInvKin libraries to calculate the inverse
% kinematics for a tensegrity spine, defined in 2 dimensions. 
% As the term is used here, 'spine'
% refers to a structure with repeated rigid bodies of the same geometry and
% mass.

% As of 2018-08-16, this scrupt (needs to) use the rigid body reformulation
% of the inverse kinematics problem. And, that reformulation currently has
% stringent constraints on when it can be used - right now, only with 2
% rigid bodies.

%% set up the workspace
clear all;
close all;
clc;

% add the core libraries, assumed to be in an adjacent folder.
addpath( genpath('../invkin_core') );
% same for the plotting.
addpath( genpath('../invkin_plotting') );

%% Set up the parameters

% debugging or not
debugging = 1;

% minimum cable force density
%q_min = 0; % units of N/m, depending on m and g
q_min = 0.5;

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 2D
% spine. Nodes are column vectors [x; z], but for ease, written as transposed
% here.

bar_endpoint = 0.5; % meters. 50 cm.

a = [   0,              0;
        bar_endpoint,   -bar_endpoint;
        -bar_endpoint,  -bar_endpoint;
        0,              bar_endpoint]';
    
if debugging
    a
end

% number of rigid bodies
b= 2;
% When removing the anchor nodes, it's like removing one of the bodies:
%b = 1;

% Configuration matrix for WHOLE STRUCTURE. (in future, pattern out.)

% Full connectivity matrix
% Rows 1-4 are cables
% Rows 5-10 are bars
% Columns 1-4 are bottom tetra nodes
% Columns 5-8 are top tetra nodes
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

if debugging
    C
end

% gravitational constant
g = 9.81;

% vector of masses of each node
% let's say each vertebra weighs 0.8 kg. Thta's about 1.7 lbs, which seems
% right to Drew if motors are included.
m_i = 0.8;

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 4, and keep nodes 5-8
% (which is vertebra two.)
%w = [0; 0; 0; 0; 1; 1; 1; 1];
% Including all nodes:
w = ones(n,1);

% In case we need it later, we can calculate the number of 'remaining'
% nodes, the not-anchored ones. Call that 'h'.
h = nnz(w);

% ...later, the command we want to create what we need is
% W = diag(w);
% W(~any(W,2), :) = [];

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
xi(4:6) = [     bar_endpoint * (3/2);
                0;
               -pi/2];
            
if debugging
    xi
end

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, z)
% calculate from position trajectory
coordinates = get_node_coordinates_2d(a, xi, debugging);

if debugging
    coordinates
end

% Reaction forces
% Hard-coded for this structure. Draw out a free body diagram to confirm
% this, or to see how you'd do it for another structure.

% get the center of mass for combined structure (don't need to split by
% rigid body.)
% Note here that I've used m_i as per-body not per-node.
% Probably accidentally changed notation w.r.t. T-CST 2018 paper.
% 4 nodes per body, so
m_node = m_i/4;

% first, the positions of each mass:
mass_positions = zeros(size(coordinates));
for i=1:n
    mass_positions(:,i) = m_node * coordinates(:,i);
end

% center of mass (com) is sum over each row, divide by total mass.
% com is 2 x 1.
com = sum(mass_positions, 2) ./ (m_node*n);

if debugging
    com
end

% Assume reaction forces in z and x at nodes 2 and 3.
% Solve AR*[R2; R3] = bR.
% AR has two components: one due to force balance, the other
% due to moment balance.
% We know the following, for our example (again, hard-coded:)
% 1) the reaction forces each act at bar_endpoint away from the origin
% (when we're summing moments around the origin)
% 2) there are no external moments in the x-direction, only z
% 3) therefore, we need a moment balance that includes all 4 reaction force
% terms (R2x, R2z, R3x, R3z), and both mg terms are around the individual
% COMs.

% ...see Drew's notes, or derive this yourself!
% abbreviation: bar_endpoint = be.
be = bar_endpoint;
% 3 rows are for sum x, sum z, moments.
AR = [1     0       1       0;
      0     1       0       1;
      be   -be     -be      -be];

bR = [0;
      2*m_i*g;
      2*m_i*g*(com(1,1))];
  
% solve. To-Do: either...
%       1) pose as an optimization problem to deal with static
%       indeterminacy when it arises
%       2) incorporate this into the optimization for the force densities?
R = AR\bR;

% R = [R2x, R2z, R3x, R3z]

% Create the vector of external forces at each node.
px = zeros(n,1);
pz = zeros(n,1);

% Only x forces are at 2 and 3.
px(2) = R(1);
px(3) = R(3);
% insert the reaction forces in z
pz(2) = R(2);
pz(3) = R(4);

% Add the gravitational reaction forces for each mass.
for i=1:n
    pz(i) = pz(i) - m_node*g;
end


%% Solve the inverse kinematics problem
%[f_opt, q_opt, Ab, pb] = invkin_core_2d_rb(x, z, px, pz, C, COMs, s, b, q_min, debugging)

% Split up the coordinates.
% The invkin core routine expects these to be column vectors, so a
% transpose is needed also
x = coordinates(1, :)';
z = coordinates(2, :)';

% For the rigid body balance, we need the COMs of each rigid body,
% not just the COM of the whole structure.
COMs = zeros(2, b);
% hard coding here: 4 nodes per body (n=8, b=2, n/b = 4.)
COMs(:,1) = sum(mass_positions(:, 1:4), 2) / (m_node*4);
COMs(:,2) = sum(mass_positions(:, 5:8), 2) / (m_node*4);

% Solve
[f_opt, q_opt, Ab, pb] = invkin_core_2d_rb(x, z, px, pz, C, COMs, s, b, q_min, w, debugging);

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
radius = 0.02; % meters.

% Plot.
plot_2d_tensegrity_invkin(C, x, z, s, radius);

















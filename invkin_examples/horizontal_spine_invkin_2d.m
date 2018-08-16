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

%% Set up the parameters

% debugging or not
debugging = 1;

% minimum cable force density
q_min = 0; % units of N/m, depending on m and g

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 2D
% spine. Nodes are column vectors [x; z], but for ease, written as transposed
% here.

bar_endpoint = 0.5 % meters. 50 cm.

a = [   0,              0;
        bar_endpoint,   -bar_endpoint;
        -bar_endpoint,  -bar_endpoint;
        0,              bar_endpoint]';
    
if debugging
    a
end

% number of rigid bodies
b= 2;
% number of nodes
n = size(a, 2) * b;

% Configuration matrix for WHOLE STRUCTURE. (in future, pattern out.)

% The nodes will be 1-5 = fixed vertebra, 6-10 = cantilevered vertebra
% The configuration sub-matrix for the bars for one vert connects the
% center node to all outer nodes
% Cv =   [1,      -1,     0,      0,      0;
%         1,      0,      -1,     0,      0;
%         1,      0,      0,      -1,     0;
%         1,      0,      0,      0,      -1];
    
% The cable connections have the vertical cables (nodes 2->7 ... 5->10)
%           1   2   3   4   5   6   7   8   9   10
% Cc_vertical =  [0,  1,  0,  0,  0,  0, -1,  0,  0,  0;
%             0,  0,  1,  0,  0,  0,  0,  -1, 0,  0;
%             0,  0,  0,  1,  0,  0,  0,  0, -1,  0;
%             0,  0,  0,  0,  1,  0,  0,  0,  0, -1];
%             
% % The saddle cables connect nodes (fixed, moving): 4+2, 4+3, 5+2, 5+3
% % adjusting for the nodes on vert 2, that's 4+7, 4+8, 5+7, 5+8
% %           1   2   3   4   5   6   7   8   9   10
% Cc_saddle =[0,  0,  0,  1,  0,  0, -1,  0,  0,  0;
%             0,  0,  0,  1,  0,  0,  0, -1,  0,  0;
%             0,  0,  0,  0,  1,  0, -1,  0,  0,  0;
%             0,  0,  0,  0,  1,  0,  0, -1,  0,  0];      

% % number of cables
% % s = 8;
% s = size(Cc_vertical, 1) + size(Cc_saddle, 1);
% % number of bars
% % r = 8;
% r = size(Cv, 1) * b;
%         
% % Finally, insert the sub-matrices block-structured into a full-size C.
% C = zeros(s + r, n);
% C(1:s, :) = [Cc_vertical; Cc_saddle];
% % rigid body connections are block-structured, so here's a shortcut with
% % kron.
% C(s+1 : end, :) = kron( eye(b), Cv);

% Full connectivity matrix
% Rows 1-4 are cables
% Rows 5-10 are bars
% Columns 1-4 are bottom tetra nodes
% Columns 5-8 are top tetra nodes
%    1  2  3  4  5  6  7  8  
C = [0  1  0  0  0 -1  0  0;  %  1
     0  0  1  0  0  0 -1  0;  %  2
     0  0  0  1  0 -1  0  0;  %  3
     0  0  0  1  0  0 -1  0;  %  4
     1 -1  0  0  0  0  0  0;  %  5
     1  0 -1  0  0  0  0  0;  %  6
     1  0  0 -1  0  0  0  0;  %  7
     0  0  0  0  1 -1  0  0;  %  8
     0  0  0  0  1  0 -1  0;  %  9
     0  0  0  0  1  0  0 -1]; % 10
 
% hard-coded here.
s = 4;
r = 6;
% ...because C is \in R^{10 x 8}.

if debugging
    C
end

% gravitational constant
g = 9.81;

% vector of masses of each node
% let's say each vertebra weighs 0.8 kg. Thta's about 1.7 lbs, which seems
% right to Drew if motors are included.
m_i = 0.8;

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
% the length of one vert is bar_endpoint. 
% x-position is coordinate 1.
xi(4:6) = [     bar_endpoint * (3/4);
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
% first, the positions of each mass:
mass_positions = zeros(size(coordinates));
for i=1:n
    mass_positions(:,i) = m_i * coordinates(:,i);
end

% center of mass (com) is sum over each row, divide by n.
% com is 2 x 1.
com = sum(mass_positions, 2) ./ n;

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
% Note here that I've used m_i as per-body not per-node.
% Probably accidentally changed notation w.r.t. T-CST 2018 paper.
% 4 nodes per body, so
m_node = m_i/4;
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

% Solve
[f_opt, q_opt, Ab, pb] = invkin_core_2d_rb(x, z, px, pz, C, com, s, b, q_min, debugging);

% For this particular problem, we get the same rank issues as with the 2D
% structure described in the Sabelhaus 2018 paper on MPC with inverse
% kinematics: i.e., that since there are 10 nodes and 3 dimensions (nd =
% 30) whereas there are only 16 members (8 cables, 4 rods per body * 2
% bodies), the A matrix is overconstrained (30 rows, 16 columns) and no
% solutions exist.

% In addition, the reduction of this problem to the rigid body model
% where instead the required A q = p has only 2bd columns (both force and
% moment for b bodies in 3 dimensions) which here is 12, but only cable
% forces are considered (s=8), still would have dimensionality issues: 
% 12 > 8, so that wouldn't necessarily work either.

% The best solution here is to allow for a set of "anchor" points where
% forces are not considered, and that only provide the coordinates for the
% force density vector, so that the fixed rigid body is not considered at
% all and its reaction forces don't need to be calculated. This would lead
% to an A q = p of (in rigid body form) 8 cables and 6 constraints (force
% moment in 3d), with the cable force densities just specified by a
% constant.

% As of 2018-07-24, such an algorithm has not been developed. It's present
% somewhat in the results from the Friesen 2014 paper although is not
% discussed at all. TO-DO: implement this.





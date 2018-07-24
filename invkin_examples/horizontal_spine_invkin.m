%% horizontal_spine_invkin.m
% Copyright Andrew P. Sabelhaus 2018

% This script used the tInvKin libraries to calculate the inverse
% kinematics for a tensegrity spine. As the term is used here, 'spine'
% refers to a structure with repeated rigid bodies of the same geometry and
% mass.

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
q_min = 2; % units of N/m, depending on m and g

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 3D
% spine. Nodes are column vectors [x; y; z], but for ease, written as transposed
% here.

bar_endpoint = 0.5 % meters. 50 cm.

a = [   0,              0,              0;
        bar_endpoint,     0,              -bar_endpoint;
        -bar_endpoint,    0,              -bar_endpoint;
        0,              bar_endpoint,     bar_endpoint;
        0,              -bar_endpoint,    bar_endpoint]';
    
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
Cv =   [1,      -1,     0,      0,      0;
        1,      0,      -1,     0,      0;
        1,      0,      0,      -1,     0;
        1,      0,      0,      0,      -1];
    
% The cable connections have the vertical cables (nodes 2->7 ... 5->10)
%           1   2   3   4   5   6   7   8   9   10
Cc_vertical =  [0,  1,  0,  0,  0,  0, -1,  0,  0,  0;
            0,  0,  1,  0,  0,  0,  0,  -1, 0,  0;
            0,  0,  0,  1,  0,  0,  0,  0, -1,  0;
            0,  0,  0,  0,  1,  0,  0,  0,  0, -1];
            
% The saddle cables connect nodes (fixed, moving): 4+2, 4+3, 5+2, 5+3
% adjusting for the nodes on vert 2, that's 4+7, 4+8, 5+7, 5+8
%           1   2   3   4   5   6   7   8   9   10
Cc_saddle =[0,  0,  0,  1,  0,  0, -1,  0,  0,  0;
            0,  0,  0,  1,  0,  0,  0, -1,  0,  0;
            0,  0,  0,  0,  1,  0, -1,  0,  0,  0;
            0,  0,  0,  0,  1,  0,  0, -1,  0,  0];      

% number of cables
% s = 8;
s = size(Cc_vertical, 1) + size(Cc_saddle, 1);
% number of bars
% r = 8;
r = size(Cv, 1) * b;
        
% Finally, insert the sub-matrices block-structured into a full-size C.
C = zeros(s + r, n);
C(1:s, :) = [Cc_vertical; Cc_saddle];
% rigid body connections are block-structured, so here's a shortcut with
% kron.
C(s+1 : end, :) = kron( eye(b), Cv);

if debugging
    C
end

% Pinned-node configuration vector: specifies which nodes are pinned joints.
% \in n x 1, where == 1 if pinned and == 0 if free.

% there are n nodes
pinned = zeros(n, 1);

% it's assumed that the fixed vertebra is built-in. Shouldn't really matter
% which of nodes 1-5 are pinned, since we don't really care about the
% reaction forces / forces inside the rigid bodies' members, so maybe we
% should check in the future if this makes a difference.
% for now, let's say all "outer" nodes are fixed.
%pinned( 2:5 ) = 1;

% four reaction forces sets was too many degrees of freedom. Instead, only
% pin nodes 2 and 3, with the center pinned now also.
% (note, the linear algebra didn't work out with just nodes 2 and 3, rank
% issues.)
pinned( 1:3 ) = 1;

% gravitational constant
g = 9.81;

% vector of masses of each node
% let's say each vertebra weighs 0.8 kg. Thta's about 1.7 lbs, which seems
% right to Drew if motors are included.
m_i = 0.8;
m = m_i * ones(n, 1);

%% Trajectory of positions

% all the positions of each rigid body (expressed as their COM positions
% and Euler angles). that's 6 states: [x; y; z; \theta, \gamma, \phi] with
% the angles being intrinsic rotations (y p r).

% use \xi for the system states.
xi = zeros(b * 6, 1);

% for rigid body 1, the fixed one, doesn't move. However, we've defined it
% in its "vertical" state, so it needs to be rotated by 90 degrees around
% the y-axis so the saddle cables can align. Let's rotate it counterclockwise so
% that nodes 4 and 5 are in +x.
xi(4:6) = [0; pi/2; 0];

% for rigid body 2, translate out in the +x direction. Translating by one
% full body length puts the tips exactly in the same plane, so maybe do it
% to 3/4 of that length.
% the length of one vert is bar_endpoint. 
% x-position is coordinate 1.
xi(7:12) = [    bar_endpoint * (3/4);
                0;
                0;
                0; % angles start here
                pi/2;
                0];
            
% Let's translate and rotate the moving vertebra 'up' a bit, to check. 
% That's 
            
if debugging
    xi
end

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, y, z)
% calculate from position trajectory
coordinates = get_node_coordinates_3d(a, xi, debugging);

if debugging
    coordinates
end

% Reaction forces
% Using the nodal positions, pinned-node configuration vector, and mass
% vector. 
% Outputs px, py, pz. *THIS ASSUMES ACTING IN GRAVITY in -Z.

[px, py, pz] = get_reaction_forces_3d(coordinates, pinned, m, g, debugging);

% slight hack for now with the sign error:
px = -px; py = -py; pz = -pz;

%% Solve the inverse kinematics problem
% invkin_core_3d(x, y, z, px, py, pz, C, s, q_min, debugging}










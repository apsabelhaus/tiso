%% horizontal_spine_invkin.m
% Copyright Andrew P. Sabelhaus 2018

% This script used the tInvKin libraries to calculate the inverse
% kinematics for a tensegrity spine. As the term is used here, 'spine'
% refers to a structure with repeated rigid bodies of the same geometry and
% mass.

%% Set up the parameters

% debugging or not
debugging = 1;

% minimum cable force density
q_min = 2; % units of N/m, depending on m and g

% Local frame for one rigid body (locations of nodes)
%a = [

% number of cables

% Configuration matrix for WHOLE STRUCTURE. (in future, pattern out.)

% Pinned-node configuration vector: specifies which nodes are pinned joints.
% \in n x 1, where == 1 if pinned and == 0 if free.

% gravitational constant
g = 9.81;

% vector of masses of each node

%% Trajectory of positions

% all the positions of each rigid body (expressed as their COM positions
% and Euler angles)

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, y, z)
% calculate from position trajectory

% Reaction forces
% Using the nodal positions, pinned-node configuration vector, and mass
% vector. 
% Outputs px, py, pz. *THIS ASSUMES ACTING IN GRAVITY.

%% Solve the inverse kinematics problem
% invkin_core_3d(x, y, z, px, py, pz, C, s, q_min, debugging}



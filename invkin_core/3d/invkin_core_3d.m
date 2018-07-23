%% invkin_core_3d.m
% Copyright Andrew P. Sabelhaus 2018

function [f_opt, q_opt, A, p] = invkin_core_3d(x, y, z, px, py, pz, C, s, q_min}
% invkin_core_3d performs a single inverse kinematics computation for a
% 3-dimensions tensegrity structure (robot.)
% This script follows the force density method for tension networks, which
% calculates a condition for static equilibrium given a specific
% structure, and calculates cable forces / force densities by minimizing
% the total force density in the structure.

% (see Schek, Tran, Friesen, etc. for formulation, and my own papers for
% more discussion on 'inverse kinematics'.)

% This function REQUIRES that the p vector already be pre-populated with
% reaction forces and gravitational forces, if appropriate.

% Inputs:
%   x, y, z = the x, y, and z positions of each node in the structure. Must
%   be the same dimension.
%
%   C = configuration matrix for the structure. See literature for
%   definition.
%
%   px, py, pz = vectors of external forces, applied to each node (has the same
%   dimension as x, y, z.)
%
%   s = number of cables (tension-only members in the system.) This is
%   required to pose the optimization program such that no cables "push."
%
%   q_min = minimum force density for an individual compression member. 

%% First, formulate the constants / parameters for the optimization:

% As aside: we can get the number of nodes and number of members by
n = size(x);
r = size(C,1) - s; % rows of C is total number of members, and s is
% provided.

% The matrix defining the length vectors between nodes
% Has dimension 
A = [ C' * diag(C * x);
      C' * diag(C * y);
      C' * diag(C * z)];
  
% Combine the p vector
p = [px; py; pz];

% Since we assume that the first s rows of C are for the cables, make the
% constraint for the min force density on those cables:
q_min_vector = q_min * [ones(s,1); zeros(r,1)];

%% Solve the optimization problem

% for now, we're going to let quadprog do the relaxation from equality
% constrained to inequality constrained.





end















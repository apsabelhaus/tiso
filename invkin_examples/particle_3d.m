%% particle_3d.m
% Copyright Andrew P. Sabelhaus 2019

% This script used the tIK library (rename later)
% to calculate the equilibrium cable tensions for the viscoelastic control
% law (Drew's dissertation work on dissipativity.)

%% set up the workspace
clear all;
close all;
clc;

% add the core libraries, assumed to be in an adjacent folder.
addpath( genpath('../invkin_core') );

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 2;

% Startup message if appropriate
if debugging >= 1
    disp('Starting particle_3D equilibrium input optimization...');
end

% minimum cable force density
q_min = 1; % units of N/m, depending on m and g

% A vector of the anchor points.
% Columns are individual points: row 1 is x, ... row 3 is z.
% Units are meters.
anchor_nodes = [0      0       -0.2    0.2;
         0.2    0.2     -0.2    -0.2;
         0.2    -0.2    0       0   ];
     
% The desired position of the point mass (the desired equil. pt.)
eq_pt = [0.05;
         0.05;
         0.05];

% Combined all points into one matrix.
nodes = [eq_pt, anchor_nodes];
    
if debugging >= 2
    nodes
end

% number of nodes
n = size(nodes, 2);
% number of cables: we have one cable per anchor, and no more.
s = size(anchor_nodes, 2);
% should be n = s + 1.

% Configuration matrix for WHOLE STRUCTURE. (in future, pattern out.)
% From dissertation: C = \begin{bmatrix} \bOnes_s & - \bI_s \end{bmatrix}
C = [ones(s,1), -eye(s)];

if debugging >= 2
    C
end

% Specifying the anchors, as points to remove from the force balance.
% Here, we only care about the particle, so
anchors_vec = [1; zeros(s,1)];

% gravitational constant
g = 9.81;

% Particle mass. The vertebra in our 2D test is about a half kilogram.
m = 0.495;

% Spring constants for each cable.
% Should be an (s \times 1) vector.
k = [300;
     100;
     150;
     350];

%% Calculations for the inputs to the core invkin library

% If we cared about reaction forces:
% [px, py, pz] = get_reaction_forces_3d(coordinates, pinned, m, g, debugging);

% ...since we really don't, just add -mg to the particle
% TO-DO: check the direction here! What convention did we use for positive
% cable forces??

px = zeros(n,1);
py = zeros(n,1);
pz = zeros(n,1);

% the particle is node 1
pz(1) = -m*g;

%% Solve the inverse kinematics problem

% Split up the coordinates.
% The invkin core routine expects these to be column vectors, so a
% transpose is needed also
x = nodes(1, :)';
y = nodes(2, :)';
z = nodes(3, :)';

% Solve
[f_opt, q_opt, A, p, lengths] = invkin_core_3d_anchorless(x, y, z, px, py, pz, C, s, q_min, anchors_vec, debugging);

% Backing out the optimal inputs,
% v_i = l_i - (F_i/kappa_i)
disp('Equilibrium inputs are:');
bar_v = lengths - (f_opt ./ k)




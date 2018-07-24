%% invkin_core_3d.m
% Copyright Andrew P. Sabelhaus 2018

function [f_opt, q_opt, A, p] = invkin_core_3d(x, y, z, px, py, pz, C, s, q_min, debugging}
% invkin_core_3d performs a single inverse kinematics computation for a
% 3-dimensional tensegrity structure (robot.)
%
% This script follows the force density method for tension networks, which
% calculates a condition for static equilibrium given a specific
% structure, and calculates cable forces / force densities by minimizing
% the total force density in the structure.
%
% IMPORTANT: this formulation includes a minimization of the forces inside
% the compression members. For rigid body robots, this might result in
% sub-optimal solutions (TBD.)
%
% (see Schek, Tran, Friesen, etc. for formulation, and my own papers for
% more discussion on 'inverse kinematics'.)
%
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
%
%   debugging = set to 1 if you want more information output on the command
%   line.

% Outputs:
%   f_opt = optimal forces (inv kin soln) for ALL members. First s are
%   cables.
%
%   q_opt = optimal force densities, corresponding to f_opt
%
%   A = calculated A matrix, provided for debugging if desired.

%% Check if inputs are valid (conformal dimensions of vectors and matrices)

% do this later...

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

if debugging
    disp('Size of A is:');
    size(A);
end

% Since we assume that the first s rows of C are for the cables, make the
% constraint for the min force density on those cables:
q_min_vector = q_min * [ones(s,1); zeros(r,1)];

%% Solve the optimization problem

% for now, we're going to let quadprog do the relaxation from equality
% constrained to inequality constrained.

if debugging
    disp('Solving the inverse kinematics problem for a 3D tensegrity structure...');
end

% quadprog's arguments are
% 0.5 * H (quadratic term), f (linear term), A_ineq, b_ineq, A_eq, b_eq

% Our problem is
% min q'q s.t. Aq=p, q - mintensionvec \geq 0

% no weighting any members different than others
% the 2 is unnecessary here since we only care about argmax
%H = 2 * eye(s + r);

% ALTERNATIVELY: just weight the cables! This way, we don't care about what
% the bars have in terms of force, so maybe there's a better minimum for
% cables-only.
% So, just have the H matrix = I for the cables, = 0 for the bars.
H = zeros(s + r);
H(1:s, 1:s) = eye(s);

% no linear term
f = [];
% no inequality constraints
A_ineq = [];
b_ineq = [];

% the equality constraints are in the form of A_eq * q <= b_eq,
% so we need q >= mintensionvec becomes -q <= -mintensionvec
A_eq = - eye(s + r);
b_eq = - q_min_vector;

% Call quadprog
[q_opt] = quadprog(H, f, A_ineq, b_ineq, A_eq, b_eq);

%% Calculate the cable forces from the force densities.

% Forces are density * length, and lengths for each member in each
% direction are
dx = C * x;
dy = C * y;
dz = C * z;

% so the lengths of each cable are the euclidean norm of each 3-vector.
% re-organize:
length_vecs = [dx'; dy'; dz'];
if debugging
    disp('Member length vectors are:');
    length_vecs
end

% the scalar lengths are then the 2-norm (euclidean) for each column, which
% is
lengths = vecnorm(length_vecs);
if debugging
    disp('Member lengths (scalar) are:');
    lengths
end

% Then, element-wise calculate the optimal forces.
f_opt = q_opt .* lengths;

end















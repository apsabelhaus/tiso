%% invkin_core_3d.m
% Copyright Andrew P. Sabelhaus 2018

function [f_opt, q_opt, A, p] = invkin_core_3d(x, y, z, px, py, pz, C, s, q_min, debugging)
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
%
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
%   debugging = Debugging level / verbosity.
%       0 = no output except for errors
%       1 = Starting message, results from quadprog
%       2 = Verbose output of status.
%
% Outputs:
%   f_opt = optimal forces (inv kin soln) for ALL members. First s are
%   cables.
%
%   q_opt = optimal force densities, corresponding to f_opt
%
%   A = calculated A matrix, provided for debugging if desired.

%% Startup message
if debugging >= 1
    disp('Starting invkin_core_3d...');
end

%% Check if inputs are valid (conformal dimensions of vectors and matrices)

% The number of nodes, according to the C matrix, is its number of columns:
n = size(C,2);
% Verify that this is also the length of all other nodal vectors.
if n ~= size(x, 1)
    error('Error: the C matrix and the x vector (node positions in x) have a different number of nodes. Cannot continue.');
elseif n ~= size(y, 1)
    error('Error: the C matrix and the z vector (node positions in z) have a different number of nodes. Cannot continue.');
elseif n ~= size(px, 1)
    error('Error: the C matrix and the px vector (node external rxn. forces in x) have a different number of nodes. Cannot continue.');
elseif n ~= size(py, 1)
    error('Error: the C matrix and the pz vector (node external rxn. forces in z) have a different number of nodes. Cannot continue.');
end

% Also check that each of these are a column vector.
if size(x, 2) ~= 1
    error('Error: x is not a column vector. Cannot continue.');
elseif size(y, 2) ~= 1
    error('Error: z is not a column vector. Cannot continue.');
elseif size(px, 2) ~= 1
    error('Error: px is not a column vector. Cannot continue.');
elseif size(py, 2) ~= 1
    error('Error: pz is not a column vector. Cannot continue.');
end

% There cannot be more cables than total number of members
if s > size(C, 1)
    error('Error: you have specified that there are more cables than total number of members, this is not possible. Cannot continue.');
end

% Minimum cable tension must be positive - otherwise we'd have cables
% holding compressive loads
if q_min < 0
    error('Error: minimum cable tension must be positive, please specify a q_min > 0.');
end

%% First, formulate the constants / parameters for the optimization:

% As aside: we can get the number of nodes and number of members by
n = size(x,1);
r = size(C,1) - s; % rows of C is total number of members, and s is
% provided.

% The matrix defining the length vectors between nodes
% Has dimension 
A = [ C' * diag(C * x);
      C' * diag(C * y);
      C' * diag(C * z)];
  
% Combine the p vector
p = [px; py; pz];

if debugging >= 2
    disp('Size of x, r, C, A is:');
    n
    r
    size(C)
    size(A)
end

% Since we assume that the first s rows of C are for the cables, make the
% constraint for the min force density on those cables:
q_min_vector = q_min * [ones(s,1); zeros(r,1)];

%% Solve the optimization problem

% for now, we're going to let quadprog do the relaxation from equality
% constrained to inequality constrained.

if debugging >= 2
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

% the inequality constraints are in the form of A_eq * q <= b_eq,
% so we need q >= mintensionvec becomes -q <= -mintensionvec
A_ineq = - eye(s + r);
b_ineq = - q_min_vector;

% equality constraints are A q = p
A_eq = A;
b_eq = p;

if debugging >= 2
    disp('Optimization parameters are:');
    H
    f
    A_eq
    b_eq
end

% Set some options for quadprog. 
% Since it's hard to understand the output in the context of our problem,
% supress quadprog's output and parse it on our own.
qp_options = optimoptions('quadprog','Display','none');

% If the user has asked for more debugging information:
if debugging >= 2
    qp_options.Display = 'final';
end

% Call quadprog
% THE BIG STEP
[q_opt, ~, exitflag, output_info] = quadprog(H, f, A_ineq, b_ineq, A_eq, ...
                                       b_eq, [], [], [], qp_options);

% Call our parser to make some sense of out what happened
if debugging >= 1
    parse_qp_invkin_output(exitflag, output_info);
end

if debugging >= 2
    disp('Optimal force densities are:');
    q_opt
end

%% Calculate the cable forces from the force densities.

% Forces are density * length, and lengths for each member in each
% direction are
dx = C * x;
dy = C * y;
dz = C * z;

% so the lengths of each cable are the euclidean norm of each 3-vector.
% re-organize:
length_vecs = [dx'; dy'; dz'];
if debugging >= 2
    disp('Member length vectors are:');
    length_vecs
end

% the scalar lengths are then the 2-norm (euclidean) for each column, which
% is
lengths = vecnorm(length_vecs);
if debugging >= 2
    disp('Member lengths (scalar) are:');
    lengths
end

% Then, element-wise calculate the optimal forces.
f_opt = q_opt .* lengths;

if debugging >= 2
    disp('Optimal forces on cables are:');
    f_opt
end

end















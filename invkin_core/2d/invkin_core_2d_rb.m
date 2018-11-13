%% invkin_core_2d_rb.m
% Copyright Andrew P. Sabelhaus 2018

function [f_opt, q_opt, Ab, pb] = invkin_core_2d_rb(x, z, px, pz, w, C, s, b, q_min, debugging)
% invkin_core_2d_rb performs a single inverse kinematics computation for a
% 2-dimensional tensegrity structure (robot), using the reformulated
% problem treating the rigid bodies with a force/moment balance and not as
% a node-graph.
%
%   *IMPORANT:* This script currently only is valid under the following
%   conditions:
%
%   d = 2 (this is a 2d script anyway!)
%   b = 2
%   Each rigid body has the same structure (specifically, same number of nodes)
%
% This script follows the force density method for tension networks, which
% calculates a condition for static equilibrium given a specific
% structure, and calculates cable forces / force densities by minimizing
% the total force density in the structure.
%
% (see Schek, Tran, Friesen, etc. for formulation, and my own papers for
% more discussion on 'inverse kinematics'.)
%
% This function REQUIRES that the p vector already be pre-populated with
% reaction forces and gravitational forces, if appropriate.

% Inputs:
%   x, z = the x, and z positions of each node in the structure. Must
%   be the same dimension.
%
%   px, pz = vectors of external forces, applied to each node (has the same
%   dimension as x, z.)
%
%   w = anchoring vector. Elements = 1 if node should be kept, = 0 if node
%   is an anchor and should be thrown out of the force balance. Size is n.
%
%   C = configuration matrix for the structure. See literature for
%   definition.
%
%   s = number of cables (tension-only members in the system.) This is
%   required to pose the optimization program such that no cables "push."
%
%   b = number of rigid bodies. With fancier algorithms, we could actually
%   do a search through the C matrix (graph search!) to find how many
%   independent cycles there are in the "r" block of it (rods only), but
%   for now, it's easier to just have the caller specify. (Robot designers
%   will know this intuitively: e.g., we've made a robot with 2 vertebrae.)
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
%   Ab, pb = calculated static equilibrium constraint, 
%           provided for debugging if desired.

%% Check if inputs are valid (conformal dimensions of vectors and matrices)

% The number of nodes, according to the C matrix, is its number of columns:
n = size(C,2);
% Verify that this is also the length of all other nodal vectors.
if n ~= size(x, 1)
    error('Error: the C matrix and the x vector (node positions in x) have a different number of nodes. Cannot continue.');
elseif n ~= size(z, 1)
    error('Error: the C matrix and the z vector (node positions in z) have a different number of nodes. Cannot continue.');
elseif n ~= size(px, 1)
    error('Error: the C matrix and the px vector (node external rxn. forces in x) have a different number of nodes. Cannot continue.');
elseif n ~= size(pz, 1)
    error('Error: the C matrix and the pz vector (node external rxn. forces in z) have a different number of nodes. Cannot continue.');
elseif n ~= size(w, 1)
    error('Error: the C matrix and the w vector (anchoring vector) have a different number of nodes. Cannot continue.');
end

% Also check that each of these are a column vector.
if size(x, 2) ~= 1
    error('Error: x is not a column vector. Cannot continue.');
elseif size(z, 2) ~= 1
    error('Error: z is not a column vector. Cannot continue.');
elseif size(px, 2) ~= 1
    error('Error: px is not a column vector. Cannot continue.');
elseif size(pz, 2) ~= 1
    error('Error: pz is not a column vector. Cannot continue.');
elseif size(w, 2) ~= 1
    error('Error: w is not a column vector. Cannot continue.');
end

% There must be at least one rigid body
if b < 1
    error('Error: b is less than one. Must be at least one rigid body.');
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
r = size(C,1) - s; % rows of C is total number of members, and s is
% provided.

% The number of nodes per rigid body is
%eta = n / b;
% Actually, now, removing the nodes that are not being considered (the
% anchors), we really want to calculate eta using h here, where h < n is
% the number of remaining nodes.
h = nnz(w);
eta = h / b;

% For later work, generalize: this is the dimensionality of the problem
% (2 dimensional.)
d = 2;

% First, let's create the A matrix.
% TO-DO: check validity for b > 2 rigid bodies.

% Remember, we're looking at something (in the end) like:
% Ab * qs = pb
% where
% Ab = [Af; Am]
% pb = [pf; pm]
% Af = forces from each cable on each rigid body
% pf = sum of external forces in each dimension for each rigid body
%       *NOTE: your calling script should include gravity here already!
% Am = sum of moments from each cable around each rigid body's center of
%       mass
% pm = sum of moments from each set of external reaction forces around each
%       rigid body's center of mass

% Helper matrix:
H_hat = [eye(s), zeros(s, r)];
% ...when multiplied by a vector of all cables, removes the bars.
% Equivalently, H_hat^\top will right-multiply the A matrices to remove the
% bars.

% NOW: We can actually collapse it in the same way as the external force
% vector. The dimensions and calculations ctually line up exactly how
% we want. Think about it like summing over chunks of columns of size eta
% within the A matrix.

% a matrix to pattern out the collapsation.
% We could do ones(1, eta) but all Drew's papers use the 'ones' as a column
% vector, so do that here, and transpose it.
collapse = kron(eye(d*b), ones(eta, 1)');

% As in the standard force-density method, calculate the static equilibrium
% matrix:
A = [ C' * diag(C * x);
      C' * diag(C * z)];
  
% ANCHORS FORMULATION:
% For using the anchors-formulation (removing the anchors from the
% constraint), we right-multiply the Af matrix.
% Each dimensional component needs to be
W = diag(w);
% Remove the zero-ed rows in W. It doesn't seem like there's a nice linear
% algebra way to do this, so use MATLAB's logical indexing instead.
W(~any(W,2), :) = [];

% This W takes out the anchors for one dimension.
% To do it in multiple dimensions, pattern it out.
Wf = kron(eye(2), W);

% We can now remove the rows from A.
% This is done "before" the rigid body reformulation,
% so this necessarily is also dependent on the assumptions: rigid bodies
% are "in a row" in terms of labeling order, etc.
Aw = Wf * A;
  
% The Af matrix can then be collapsed to size (2b x s) from size (2n x
% (s+r))
Af = collapse * Aw * H_hat';

if debugging
    A
    Aw
    Af
end
  
% Combine the p vector
p = [px; pz];

% We need to remove the not-used nodes from the external balance now too.
% The Wf matrix works the same way here (nicely same dimensions.)
pw = Wf * p;

% Now, we can also calculate pf by "collapsing" it down to per-body instead
% of per-node.
pf = collapse * pw;

% NOTE: As of Nov. 12th, have substituted the general moment balance
% formulation. This seems to be correct, much more so than the old
% iterative calculation (which was always suspect.) 

% The moment balance: 
% We can sum around the origin, and then remove rows similarly to Aw.
% For example, for the moments due to external forces, we do something
% like:
% M_ext = det([x, y; Fx, Fy])
%       = x*Fy - y*Fx
% where x and y are the node positions at which [Fx, Fy] act.
% This is equivalent to multiplying the p vector by the matrix
% TO-DO interchange z with y since right-hand convention means y is "up" in
% 2D.
B = [-diag(z), diag(x)];

% the B matrix is \in \mathbb{R}^(n, 2n} since we're only doing moment
% balance in the Z direction.

% The moments from the total external force at each node are
pm = B * p;

% Then collapse down according to (a) removing anchors and (b) rigid body
% reformulation.

% What's the 'collapse' matrix to use here? It's different dimensions.
% Need to sum according to eta, but x and y are combined, so should be
% half the number of columns (only n nodes, not 2n) and half the number of
% rows (only summing in Z for each body, not summing and X and Y for each
% body)
% It would then be
collapse_m = kron(eye(b), ones(eta, 1)');

% The moment contribution from external forces should be size b x 1,
% one sum moment in Z for each rigid body.
% We also need to remove the anchors, as with the force balance.
pm = collapse_m * W * pm;

% For the collapsation of the forces on the lefthand side, due to the
% cables, we perform exactly the same operation. 
% Multiplying by B calculates the moments due to each cable around the
% origin of the world frame
Am = B * A;

% Then remove the bars, and collapse it down to one moment balance.
% This leaves a matrix of size b x s.
Am = collapse_m * W * Am * H_hat';

% finally, can concatenate the force balance with the moment balance.
Ab = [Af; Am];
pb = [pf; pm];

if debugging
    disp('Size of x, r, C, Ab, pb is:');
    n
    r
    size(C)
    size(Ab)
    size(pb)
    disp('Values of Ab and pb are:');
    Ab
    pb
end

% Since we assume that the first s rows of C are for the cables, make the
% constraint for the min force density on those cables:
q_min_vector = q_min * [ones(s,1)];

%% Solve the optimization problem

% for now, we're going to let quadprog do the relaxation from equality
% constrained to inequality constrained.

if debugging
    disp('Solving the inverse kinematics problem for a 2D tensegrity structure...');
end

% quadprog's arguments are
% 0.5 * H (quadratic term), f (linear term), A_ineq, b_ineq, A_eq, b_eq

% Our problem is
% min qs'qs s.t. Ab*qs=pb, q - mintensionvec \geq 0

% no weighting any members different than others
% the 2 is unnecessary here since we only care about argmax
%H = 2 * eye(s + r);

% Our H is only in s x s now
H = eye(s);

% no linear term
f = [];

% the inequality constraints are in the form of A_eq * q <= b_eq,
% so we need q >= mintensionvec becomes -q <= -mintensionvec
A_ineq = - eye(s);
b_ineq = - q_min_vector;

% equality constraints are A q = p
A_eq = Ab;
b_eq = pb;

if debugging
    disp('Optimization parameters are:');
    H
    f
    A_eq
    b_eq
end

% TO-DO: sanitize here? Maybe zero-out anything that hits numberical
% precision, e.g. any element < epsilon, where epsilon = 1e-14 or so,
% set = 0. Is this a good idea?

% Call quadprog
[q_opt] = quadprog(H, f, A_ineq, b_ineq, A_eq, b_eq);

if debugging
    disp('Optimal force densities are:');
    q_opt
end

%% Calculate the cable forces from the force densities.

% % Forces are density * length, and lengths for each member in each
% % direction are
% dx = C * x;
% dz = C * z;

% all the lengths of each cable are:
dx = H_hat * C * x;
dz = H_hat * C * z;

% so the lengths of each cable are the euclidean norm of each 2-vector.
% re-organize:
length_vecs = [dx'; dz'];
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

% Then, calculate the optimal forces. (results should be same dim as
% q_opt.)
% again, probably some nice linear algebra here (TO-DO: efficiency),
% but looping is intuitive and I'm busy.
f_opt = zeros(size(q_opt));
for k=1:s
    f_opt(k) = lengths(k) * q_opt(k);
end
%f_opt = lengths * q_opt;

if debugging
    disp('Optimal forces on cables are:');
    f_opt
end

end















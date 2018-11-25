%% invkin_core_2d_rb.m
% Copyright Andrew P. Sabelhaus 2018

function [f_opt, q_opt, Ab, pb] = invkin_core_2d_rb(x, y, px, py, w, C, s, b, q_min, debugging)
% invkin_core_2d_rb performs a single inverse kinematics computation for a
% 2-dimensional tensegrity structure (robot), using the reformulated
% problem treating the rigid bodies with a force/moment balance and not as
% a node-graph.
%
%   *IMPORANT:* This script currently only is valid under the following
%   conditions:
%
%   d = 2 (this is a 2d script anyway!)
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
% However, if using the anchor-removal formulation, and all nodes with
% external reaction forces are removed via w, then this no longer matters
% and a zero vector can be passed in.
%
% Inputs:
%   x, y = the x, and y positions of each node in the structure. Must
%   be the same dimension. Cartesian frame.
%       Note, in two dimensions, y is "up" e.g. gravity works in -y.
%       This is so that a right-handed coordinate frame still makes sense
%       in 2D.
%
%   px, py = vectors of external forces, applied to each node (has the same
%   dimension as x, y.)
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
%   debugging = level of debugging/verbosity from this script. Options:
%       0 = no output except for errors
%       1 = Starting message, results from quadprog
%       2 = Verbose output, lots of dimensions and matrices and whatnot.
%
% Outputs:
%   f_opt = optimal forces (inv kin soln) for ALL members. First s are
%   cables.
%
%   q_opt = optimal force densities, corresponding to f_opt
%
%   Ab, pb = calculated static equilibrium constraint, 
%           provided for debugging if desired.

%% Startup message.

if debugging >= 1
    disp('Starting invkin_core_2d_rb...');
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
elseif n ~= size(w, 1)
    error('Error: the C matrix and the w vector (anchoring vector) have a different number of nodes. Cannot continue.');
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
% anchors), we really want to calculate eta using h here, where h <= n is
% the number of remaining nodes.
% THIS WORKS EVEN IF w = ones, since in that case, h = n.
h = nnz(w);
eta = h / b;

% A warning message here, since it's easy to mess up with the anchor
% removal code, and since no dimensional errors are generated if so.
if debugging >= 1
    if h < n
        disp('You have removed anchored nodes from the formulation.');
        disp('Remember to change b, the number of bodies, or else you will get errors.');
    end
end

% For later work, generalize: this is the dimensionality of the problem
% (2 dimensional.)
d = 2;

% First, let's create the A matrix.

% Remember, we're looking at something (in the end) like:
% Ab * qs = pb
% where
% Ab = [Af; Am]
% pb = [pf; pm]
% Af = forces from each cable on each rigid body
% pf = sum of external forces in each dimension for each rigid body
%       *NOTE: your calling script should include gravity here already!
% Am = sum of moments from each cable (we sum around origin, it's easier)
% pm = sum of moments from each set of external reaction forces around the
%       origin (easier than CoM and still valid.)

% Helper matrix:
H_hat = [eye(s), zeros(s, r)];
% ...when multiplied by a vector of all cables, removes the bars.
% Equivalently, H_hat^\top will right-multiply the A matrices to remove the
% bars.


% A matrix to pattern out the collapsation.
% We could do ones(1, eta) but all Drew's papers use the 'ones' as a column
% vector, so do that here, and transpose it.
K = kron(eye(d*b), ones(eta, 1)');

% As in the standard force-density method, calculate the static equilibrium
% matrix:
A = [ C' * diag(C * x);
      C' * diag(C * y)];
  
% ANCHORS FORMULATION:
% For using the anchors-formulation (removing the anchors from the
% constraint), we left-multiply the Af matrix.
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
Af = K * Aw * H_hat'; 
% Af = K * A * H_hat';

% We can do a quick check: if any individual column is zero, that means
% that a cable somehow has zero contribution to the force balance, which is
% not possible and means that the user specified some incorrect combination
% of w and b.
% Quick MATLAB trick:
if any(~any(Af))
    error('Error, some columns in Af are zero. This likely means that you have removed anchors without changing b, the number of bodies, or otherwise specified incorrect parameters. Exiting.');
end

if debugging >= 2
    disp('Equilibrium constraint, force:');
    A
    Aw
    Af
end
  
% Combine the p vector
p = [px; py];

% We need to remove the not-used nodes from the external balance now too.
% The Wf matrix works the same way here (nicely same dimensions.)
pw = Wf * p;

% Now, we can also calculate pf by "collapsing" it down to per-body instead
% of per-node.
pf = K * pw;

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
% This is equivalent to multiplying the p vector by the matrix:
B = [-diag(y), diag(x)];

% the B matrix is \in \mathbb{R}^(n, 2n} since we're only doing moment
% balance in the Z direction.

% The moments from the total external force at each node are:
% (this vector is size n, since we only have one total moment at each node,
% not three as in the 3D case.)
p_moments = B * p;

% Then collapse down according to (a) removing anchors and (b) rigid body
% reformulation.

% What's the 'collapse' matrix to use here? It's different dimensions.
% Need to sum according to eta, but x and y are combined, so should be
% half the number of columns (only n nodes, not 2n) and half the number of
% rows (only summing in Z for each body, not summing and X and Y for each
% body)
% It would then be
K_m = kron(eye(b), ones(eta, 1)');

% The moment contribution from external forces should be size b x 1,
% one sum moment in Z for each rigid body.
% We also need to remove the anchors, as with the force balance.
pm = K_m * W * p_moments;

% For the collapsation of the forces on the lefthand side, due to the
% cables, we perform exactly the same operation. 
% Multiplying by B calculates the moments due to each cable around the
% origin of the world frame
Am = B * A;

% Then remove the bars, and collapse it down to one moment balance.
% This leaves a matrix of size b x s.
Am = K_m * W * Am * H_hat';

% finally, can concatenate the force balance with the moment balance.
Ab = [Af; Am];
pb = [pf; pm];

if debugging >= 2
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

if debugging >= 2
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

if debugging >= 2
    disp('Optimization parameters are:');
    H
    f
    A_eq
    b_eq
end

% TO-DO: sanitize here? Maybe zero-out anything that hits numberical
% precision, e.g. any element < epsilon, where epsilon = 1e-14 or so,
% set = 0. Is this a good idea?

% Set some options for quadprog. 
% Since it's hard to understand the output in the context of our problem,
% supress quadprog's output and parse it on our own.
qp_options = optimoptions('quadprog','Display','none');

% If the user has asked for more debugging information:
if debugging >= 2
    qp_options.Display = 'final';
end

% Call quadprog
% THE BIG STEP!
[q_opt, ~, exitflag, output_info] = quadprog(H, f, A_ineq, b_ineq, A_eq, ...
                                       b_eq, [], [], [], qp_options);

% Call our parser to make some sense of out what happened
if debugging >= 1
    parse_qp_invkin_output(exitflag, output_info);
end

% more output if higher level of verbosity
if debugging >= 2
    disp('Optimal force densities are:');
    q_opt
end

% We also need to check, afterwards, if the constraints were satisfied.
% In particular, the 'max iterations exceeded' may emit solutions that are
% OK-enough but also may emit invalid solutions.
% So, quick check, independent of debugging level.
% If any element of Ab * q_opt does not equal pb,
% e.g. they are outside some tolerance level epsilon,
eps = 1e-10;
if any( (Ab * q_opt - pb) > eps )
    warning(strcat('invkin_core_2d_rb returned a result but the constraints were NOT satisfied to a tolerance of ', ...
        num2str(eps), ', solution may be invalid.'));
end


%% Calculate the cable forces from the force densities.

% % Forces are density * length, and lengths for each member in each
% % direction are
% dx = C * x;
% dz = C * z;

% all the lengths of each cable are:
dx = H_hat * C * x;
dz = H_hat * C * y;

% so the lengths of each cable are the euclidean norm of each 2-vector.
% re-organize:
length_vecs = [dx'; dz'];
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

% Then, calculate the optimal forces. (results should be same dim as
% q_opt.)
% again, probably some nice linear algebra here (TO-DO: efficiency),
% but looping is intuitive and I'm busy.
f_opt = zeros(size(q_opt));
for k=1:s
    f_opt(k) = lengths(k) * q_opt(k);
end
%f_opt = lengths * q_opt;

if debugging >= 2
    disp('Optimal forces on cables are:');
    f_opt
end

end















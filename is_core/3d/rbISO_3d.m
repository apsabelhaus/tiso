%% rbISO_3d.m
% Copyright Andrew P. Sabelhaus and Albert Li 2019

function [fOpt, qOpt, lengths, Ab, pb] = rbISO_3d(mandInputs, varargin)
%% rbISO_3d  
% rbISO_3d performs an inverse statics optimziation for a single pose of a
% 3-dimensional tensegrity structure (robot), using the reformulated
% rigid body (rb) force/moment balance.
%
%   *IMPORANT:* This script currently only is valid under the following
%   conditions:
%
%   Each rigid body has the same number of nodes
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
%
% Inputs:
%
%   mandInputs: a struct of inputs to the problem that are MANDATORY. These
%       specifically include:
%
%           x, y, z = the x, y, and z positions of each node in the structure. Must
%               be the same dimension. Cartesian frame.
%               Note, in three dimensions, z is "up" e.g. gravity works in -z.
%
%           px, py, pz = vectors of external forces, applied to each node (has the same
%               dimension as x, y, z.)
%
%           C = configuration matrix for the structure. See literature for
%               definition.
%
%           s = number of cables (tension-only members in the system.) This is
%               required to pose the optimization program such that no cables "push."
%
%           b = number of rigid bodies. If eta is not specified as an
%               optional argument, then it's assumed that there are the
%               same number of nodes in all bodies, eta = n/b (or adjusted
%               by W, to h/b.)
%
%   optionalInputs: a struct of inputs to the problem that are OPTIONAL.
%       These specifically include (and have the default values/behaviors
%       of):
%
%           w = anchoring vector. Elements = 1 if node should be kept, = 0 if node
%               is an anchor and should be thrown out of the force balance. Size is n.
%               Default value: ones(n,1), keeping all nodes.
%
%           R = weighting matrix for the optimization problem. Use helper functions
%               to choose it (either 2-norm of force densities, or if spring constants
%               available, total potential energy in the cables.)
%               Default value: eye(s), i.e., getObj_2norm.
%
%           qMin = minimum force density for an individual tension member. 
%               Scalar. (TO-DO: implement qMin possible size s.)
%               Default value: 0.
%
%           debugging = level of debugging/verbosity from this function. Options:
%               0 = no output except for errors
%               1 = Starting message, results from quadprog
%               2 = Verbose output, lots of dimensions and matrices and whatnot.
%               Default value: 1.
%
%           kappa = spring constant(s) for the cables. Function is polymorphic
%               with either:
%                   \in \mathbb{R}_+ : one spring constant for all cables
%                   \in \mathbb{R}_+^s : vector of different spring constants
%                                        per-cable
%               Default value: unused, only needed for minimum rest length
%               constraint (which is optional).
%
%           uMin = minimum cable rest length for an individual cable. Scalar, \in
%               \mathbb{R}^+.
%               Default value: unused, only needed for minimum rest length
%               constraint (which is optional.)
%
%           eta \in \mathbb{Z}_+^(b) = a vector of number of nodes in a body, 
%               in each entry is the number of nodes for that body. This
%               allows specifying unevenly-distributed numbers of nodes
%               (not just n/b.)
%
% Outputs:
%   fOpt = optimal forces (inv kin soln) for ALL members. First s are
%   cables.
%
%   qOpt = optimal force densities, corresponding to fOpt
%
%   lengths = lengths of each cable, as specified via C and {x,y,z}. This
%       could be calc'd elsewhere since it doesn't actually depend on the
%       output of quadprog, but it's useful to return for calculating u
%       from q* in the caller.
%
%   Ab, pb = calculated static equilibrium constraint, 
%           provided for debugging if desired.

%% Check if inputs are valid (conformal dimensions of vectors and matrices)

% First: do we have optional arguments/inputs?
hasOI = (nargin > 1);

% pick out the optional args struct if it's there
if hasOI
    optionalInputs = varargin{1};
else
    % declaring a dummy here helps with assigning default values later
    optionalInputs = struct();
end

% validate each way. This is so that we don't pass to validate_2d a cell
% array of size one, containing one empty cell array.
if hasOI
    validate_3d(mandInputs, optionalInputs);
else
    validate_3d(mandInputs);
end

% Assuming validation works, expand out mandatory arguments
x = mandInputs.x;
y = mandInputs.y;
z = mandInputs.z;
px = mandInputs.px;
py = mandInputs.py;
pz = mandInputs.pz;
C = mandInputs.C;
s = mandInputs.s;
b = mandInputs.b;
% The number of nodes, according to the C matrix, is its number of columns:
n = size(C,2);

% get the default optional arguments to use / compare against those passed
% in. A struct.
defaultOI = getDefaultOptInputs(n, s);

% expand optional arguments if they exist, fill in defaults for those that
% do not

if isfield(optionalInputs, 'w')
    w = optionalInputs.w;
else
    w = defaultOI.w;
end

if isfield(optionalInputs, 'R')
    R = optionalInputs.R;
else
    R = defaultOI.R;
end
   
if isfield(optionalInputs, 'qMin')
    qMin = optionalInputs.qMin;
else
    qMin = defaultOI.qMin;
end
    
if isfield(optionalInputs, 'debugging')
    debugging = optionalInputs.debugging;
else
    debugging = defaultOI.debugging;
end
    
if isfield(optionalInputs, 'kappa')
    kappa = optionalInputs.kappa;
    % Expand kappa if a scalar was passed in.
    if isscalar(kappa)
        kappa = kappa * ones(s,1);
    end
    % else do nothing - spring const not defined
end

if isfield(optionalInputs, 'uMin')
    uMin = optionalInputs.uMin;
    % else do nothing - spring const not defined
end

if isfield(optionalInputs, 'eta')
    eta = optionalInputs.eta;
end

% default is to not use the minimum rest length constraint
uMinConstr = 0;
% but if both uMin and kappa are defined, use it:
if (exist('kappa', 'var') == 1) && (exist('uMin', 'var') == 1)
    uMinConstr = 1;
end


%% Startup message.

if debugging >= 1
    disp('Starting rbISO_3d (inverse statics optimization, three dimensions, rigid body reformulation)...');
end

%% First, formulate the constants / parameters for the optimization:

% As aside: we can get the number of nodes and number of members by
r = size(C,1) - s; % rows of C is total number of members, and s is
% provided.

% For below, we're going to now (as of 2019-07-24) adopt the notation that
% eta is a column vector of size 'b', with number of nodes in each body in
% it.

% Previous behavior is recovered here, if eta is not specified, by assuming
% that the nodes are evenly distributed in the remaining bodies.

% We really want to calculate eta using h here, where h <= n is
% the number of remaining nodes.
% THIS WORKS EVEN IF w = ones, since in that case, h = n.
h = nnz(w);

if ~exist('eta', 'var')
    nodes_per_body = h / b;
    eta = ones(b,1) * nodes_per_body;
end

% A warning message here, since it's easy to mess up with the anchor
% removal code, and since no dimensional errors are generated if so.
if debugging >= 1
    if h < n
        disp('You have removed anchored nodes from the formulation.');
        disp('Remember to change b, the number of bodies, or else you will get errors.');
    end
end

% For later work, generalize: this is the dimensionality of the problem
% (3 dimensional.)
d = 3;

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
H = get_H(s, r);
 
% ...when multiplied by a vector of all cables, removes the bars:
% q_s = H^\top q.
% Equivalently, H will right-multiply the A matrices to remove the bars.

% For later, calculate the lengths.
lengths = getLen_3d(x, y, z, s, C);

% A matrix to pattern out the collapsation.
% We could do ones(1, eta) but all Drew's papers use the 'ones' as a column
% vector, so do that here, and transpose it.
% K = kron(eye(d*b), ones(eta, 1)');

% Use the new helper function.
[~, K] = get_Kb(eta, d);

% As in the standard force-density method, calculate the static equilibrium
% matrix:
A = [ C' * diag(C * x);
      C' * diag(C * y);
      C' * diag(C * z)];
  
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
Wf = kron(eye(d), W);

% We can now remove the rows from A.
% This is done "before" the rigid body reformulation,
% so this necessarily is also dependent on the assumptions: rigid bodies
% are "in a row" in terms of labeling order, etc.
Aw = Wf * A;
  
% The Af matrix can then be collapsed to size (db x s) from size (dn x
% (s+r))
Af = K * Aw * H; 

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
p = [px; py; pz];

% We need to remove the not-used nodes from the external balance now too.
% The Wf matrix works the same way here (nicely same dimensions.)
pw = Wf * p;

% Now, we can also calculate pf by "collapsing" it down to per-body instead
% of per-node.
pf = K * pw;

% The moment balance: 
% We can sum around the origin, and then remove rows similarly to Aw.
% Using the helper function,
B = getB_3d(n, x, y, z);

% The moments from the total external force at each node are:
% (this vector is size n, since we only have one total moment at each node,
% not three as in the 3D case.)
p_moments = B * p;

% Then collapse down according to (a) removing anchors and (b) rigid body
% reformulation.

% Using the same K as above,
pm = K * Wf * p_moments;

% For the collapsation of the forces on the lefthand side, due to the
% cables, we perform exactly the same operation. 
% Multiplying by B calculates the moments due to each cable around the
% origin of the world frame
Am = B * A;

% Then remove the bars, and collapse it down to one moment balance.
Am = K * Wf * Am * H;

% finally, can concatenate the force balance with the moment balance.
Ab = [Af; Am];
pb = [pf; pm];

if debugging >= 2
    disp('n, r, C, Ab, pb is:');
    n
    r
    size(C)
    size(Ab)
    size(pb)
    disp('Values of Ab and pb are:');
    Ab
    pb
end

% the inequality constraints are in the form of A_ineq * q <= b_ineq,
% so we need q >= mintensionvec becomes -q <= -mintensionvec
% Since we assume that the first s rows of C are for the cables, make the
% constraint for the min force density on those cables:
qMinVec = qMin * ones(s,1);

% If including the minimum rest length constraint, this becomes

if uMinConstr == 1
    % The inequality constraint for the input saturation constraint is
    % \mathbf{L} \mathbf{q} \leq \bm{\kappa}(\boldsymbol{\ell} -
    % \mathbf{u}^{min}) so
    Kappa = diag(kappa);
    L = diag(lengths);
    uMinVec = uMin * ones(s,1);
    % so the righthandside becomes
    bIneqSat = Kappa * (lengths - uMinVec);
    
    % the total inequality constraint is then
    A_ineq = [L;
             -eye(s)];
    b_ineq = [ bIneqSat;
              -qMinVec];
else
    % If not using the input saturation constraint,
    A_ineq = - eye(s);
    b_ineq = - qMinVec;
end


%% Solve the optimization problem

% for now, we're going to let quadprog do the relaxation from equality
% constrained to inequality constrained.

if debugging >= 2
    disp('Solving the inverse statics problem for a 3D tensegrity structure...');
end

% quadprog's arguments are
% 0.5 * R (quadratic term), f (linear term), A_ineq, b_ineq, A_eq, b_eq

% no linear term
f = [];

% equality constraints are A q = p
A_eq = Ab;
b_eq = pb;

if debugging >= 2
    disp('Optimization parameters are:');
    A_ineq  
    b_ineq
    A_eq
    b_eq
end

% TO-DO: sanitize here? Maybe zero-out anything that hits numberical
% precision, e.g. any element < epsilon, where epsilon = 1e-14 or so,
% set = 0. Is this a good idea?

% Set some options for quadprog. 
% Since it's hard to understand the output in the context of our problem,
% supress quadprog's output and parse it on our own.
qpOptions = optimoptions('quadprog','Display','none');

% If the user has asked for more debugging information:
if debugging >= 2
    qpOptions.Display = 'final';
end

% Call quadprog
% THE BIG STEP!
[qOpt, ~, exitFlag, outputInfo] = quadprog(R, f, A_ineq, b_ineq, A_eq, ...
                                       b_eq, [], [], [], qpOptions);

% Call our parser to make some sense of out what happened
if debugging >= 1
    parseISOresult(exitFlag, outputInfo);
end

% more output if higher level of verbosity
if debugging >= 2
    disp('Optimal force densities are:');
    qOpt
end

% We also need to check, afterwards, if the constraints were satisfied.
% In particular, the 'max iterations exceeded' may emit solutions that are
% OK-enough but also may emit invalid solutions.
% So, quick check, independent of debugging level.
% If any element of Ab * qOpt does not equal pb,
% e.g. they are outside some tolerance level epsilon,
eps = 1e-10;
if any( (Ab * qOpt - pb) > eps )
    warning(strcat('rbISO_2d returned a result but the constraints were NOT satisfied to a tolerance of ', ...
        num2str(eps), ', solution may be invalid.'));
end


%% Calculate the cable forces from the force densities.

% Then, calculate the optimal forces. (results should be same dim as
% qOpt.)
% again, probably some nice linear algebra here (TO-DO: efficiency),
% but looping is intuitive and I'm busy.
fOpt = zeros(size(qOpt));
for k=1:s
    fOpt(k) = lengths(k) * qOpt(k);
end

if debugging >= 2
    disp('Optimal forces on cables are:');
    fOpt
end

end















% Andrew P. Sabelhaus and Albert H. Li, 2019

function validate_3d(mi, varargin)
%% validate_3d
%
%   Validates the inputs passed into an ISO problem.
%   ****As of June 2019: works for rigid body reformulation ONLY.
   
%% extracting mandatory inputs:
   
%           x, y, z = the x, and y, z positions of each node in the structure. Must
%               be the same dimension. Cartesian frame.
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
%           b = number of rigid bodies. With fancier algorithms, we could actually
%               do a search through the C matrix (graph search!) to find how many
%               independent cycles there are in the "r" block of it (rods only), but
%               for now, it's easier to just have the caller specify. (Robot designers
%               will know this intuitively: e.g., we've made a robot with 2 vertebrae.)
%

if(any(structfun(@isempty,mi)) == 1)
    % check if a field is defined but doesn't have a value.
   error('You are missing mandatory inputs!')
end

x = mi.x;
y = mi.y;
z = mi.z;
px = mi.px;
py = mi.py;
pz = mi.pz;
C = mi.C;
s = mi.s;
b = mi.b;
    
%% Validate the mandatory arguments

%%%%%% TO-DO: verify integers vs. doubles, scalar, ...

% The number of nodes, according to the C matrix, is its number of columns:
n = size(C,2);

% s and b should be scalars
if ~isscalar(s)
    error('Error: s is not scalar.');
    
elseif ~isscalar(b)
    error('Error: b is not scalar.');
    
end

% Verify that n is also the length of all other nodal vectors.
if n ~= size(x, 1)
    error('Error: the C matrix and the x vector (node positions in x) have a different number of nodes. Cannot continue.');
    
elseif n ~= size(y, 1)
    error('Error: the C matrix and the y vector (node positions in y) have a different number of nodes. Cannot continue.');
    
elseif n ~= size(z, 1)
    error('Error: the C matrix and the z vector (node positions in z) have a different number of nodes. Cannot continue.');
    
elseif n ~= size(px, 1)
    error('Error: the C matrix and the px vector (node external rxn. forces in x) have a different number of nodes. Cannot continue.');
    
elseif n ~= size(py, 1)
    error('Error: the C matrix and the py vector (node external rxn. forces in y) have a different number of nodes. Cannot continue.');
    
elseif n ~= size(pz, 1)
    error('Error: the C matrix and the pz vector (node external rxn. forces in z) have a different number of nodes. Cannot continue.');

% Also check that each of these are a column vector.
elseif size(x, 2) ~= 1
    error('Error: x is not a column vector. Cannot continue.');
    
elseif size(y, 2) ~= 1
    error('Error: y is not a column vector. Cannot continue.');
    
elseif size(z, 2) ~= 1
    error('Error: z is not a column vector. Cannot continue.');
    
elseif size(px, 2) ~= 1
    error('Error: px is not a column vector. Cannot continue.');
    
elseif size(py, 2) ~= 1
    error('Error: py is not a column vector. Cannot continue.');
    
elseif size(pz, 2) ~= 1
    error('Error: pz is not a column vector. Cannot continue.');
    
end

% There must be at least one rigid body
if b < 1
    error('Error: b is less than one. Must be at least one rigid body.');
end

% There cannot be more cables than total number of members
if s > size(C, 1)
    error('Error: you have specified that there are more cables than total number of members, this is not possible. Cannot continue.');
end

%% Validate the optional arguments

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

% First: do we have optional arguments/inputs?
hasOI = (nargin > 1);

% pick out the optional args struct if it's there
if hasOI
    oi = varargin{1};
else
    oi = struct();
end


% In order, pick out each optional argument, and do a check on it.
if isfield(oi, 'w')
    w = oi.w;
    if size(w, 2) ~= 1
        error('Error: w is not a column vector. Cannot continue.');
    elseif n ~= size(w, 1)
        error('Error: the C matrix and the w vector (anchoring vector) have a different number of nodes. Cannot continue.');
    end

elseif isfield(oi, 'R')
    R = oi.R;
    % to-do: fill in
   
elseif isfield(oi, 'qMin')
    qMin = oi.qMin;
    % Minimum cable tension must be positive - otherwise we'd have cables
    % holding compressive loads
    if qMin < 0
        error('Error: minimum cable tension must be positive, please specify a qMin > 0.');
    end
    
elseif isfield(oi, 'debugging')
    debugging = oi.debugging;
    % to-do: fill in
    
elseif isfield(oi, 'kappa')
    kappa = oi.kappa;
    % validation for the spring constants
    if ~isscalar(kappa)
        % then it needs to have one entry per cable
        if size(kappa, 1) ~= s
            error('Error: spring constant vector kappa is not correctly sized. Must be scalar (all cables same) or vector of size s.');
        end
    elseif size(kappa, 2) ~= 1
        error('Error: spring constant vector is not a column vector.');
    elseif any(kappa <= 0)
        error('Error: nonpositive spring constants! Must be > 0.');
    end

elseif isfield(oi, 'uMin')
    uMin = oi.uMin;
    % to-do: fill in
    
end

   
end

% getFrictionlessRxn_3d.m
% Copyright Andrew P. Sabelhaus 2019

function [rx, ry, rz] = getFrictionlessRxn_3d(x, y, z, pinned, m, g, debugging)
% getFrictionlessRxn_3d.m
%   Calculates reaction forces at pinned joints in 3d with the constraint
%   that rx = 0 and ry = 0, only allowing forces in Z. I.e., wanting a
%   reaction force balance with no friction allowed at the pinned joints.
%
%   Inputs:
%       x, y, z \in \mathbb{R}^n, coordinates of all nodes
%
%       pinned \in {0,1}^n, vector of 0s and 1s designating which nodes
%       should have reaction forces in z
%
%       m \in \mathbb{R}^n, mass allocated at each node
%
%       g, grav constant
%
%       debugging, debugging level parameter (0, 1, 2)
%
%   Outputs:
%
%       rx, ry, rz \in \mathbb{R}^n, reaction forces due to global
%       force/moment balance. These should be 0 for rx and ry and nonzero
%       in rz ONLY for nodes that have pinned = 1 at that node.
%

% Setup: pull out some constants
% total number of nodes
n = size(x, 1);
% total number of pinned nodes
v = nnz( pinned == 1 );

if( debugging >= 2)
    disp('Number of nodes and pinned nodes:');
    n
    v
end

% Pick out the x, y, z coords for the pinned joints only.
% start with all the coordinates,
xv = x;
yv = y;
zv = z;
% Remove the unpinned rows (this comes from the W calculation in the rbISO
% scripts)
xv( ~any(pinned,2) ) = [];
yv( ~any(pinned,2) ) = [];
zv( ~any(pinned,2) ) = [];

if( debugging >= 2)
    disp('Nodal vectors for pinned joints:');
    xv
    yv
    zv
end

% The constraint is made up of three blocks:
% (1) frictionless constraint (all rx, ry set to zero)
% (2) force balance constraint
% (3) moment balance constraint

% (1) frictionless constraint
% no easy way to 'kron' this, so write it out explicitly
A_fric = [eye(v),   zeros(v),   zeros(v);
          zeros(v), eye(v),     zeros(v)];
% no frictional forces allowed
b_fric = zeros(2*v, 1);

if( debugging >= 2)
    disp('Constraint for frictionless is:');
    A_fric
    b_fric
end

% (2) force balance constraint
% we also don't need node locations here, just sum over the whole structure
% collapsing for all pinned nodes, one constraint per direction, whole
% structure
Kv = kron(eye(3), ones(v,1)');
A_fbal = Kv;
% A_fbal = kron(eye(3), ones(v,1)');
% reaction forces counteract gravity in the Z direction
% NOTE that this is NOT -mg, since it's Rxn - mg = 0 becomes Rxn = mg
% abuse MATLAB's broadcasting, m \in \mathbb{R}^n already
grav_forces = m * g;
% external forces on all nodes due to grav
p_ext = [zeros(n,1); zeros(n,1); grav_forces];
% combine in all 3 directions (x, y, z)
% this is for ALL NODES not the pinned ones (all nodes have grav)
Kn = kron(eye(3), ones(n,1)');
% combined 
b_fbal = Kn * p_ext;

if( debugging  >= 2 )
    disp('Constraint for force balance is:');
    A_fbal
    b_fbal
end

% (3) moment balance constraint
% this is where we finally use the nodal coordinate vectors
% this follows very closely to the rbISO moment balance

% moment arms for the pinned nodes
% we've already got a function for this cross product matrix
Bv = getB_3d(v, xv, yv, zv);
% collapsing moments for whole structure in all directions
A_mbal = Kv * Bv;

% moment arms for ALL nodes (all grav forces)
Bn = getB_3d(n, x, y, z);
% collapsing moments for whole structure in all directions
b_mbal = Kn * Bn * p_ext;

if( debugging >= 2 )
    disp('Constraint for moment balance is:');
    A_mbal
    b_mbal
end

% Finally, stack everything into one system of equations,
A = [A_fric; A_fbal; A_mbal];
b = [b_fric; b_fbal; b_mbal];

% the pseudoinverse gives the least squares solution, but need to check if
% the system even has solutions first
existence = ( rank([A, b]) >= rank(A) );

if( debugging >= 2 )
    disp('Solutions to getFrictionlessRxn_3d exist?')
    existence
    disp('Compare rank test to number unknowns, if greater than 0, infinitely many solutions:')
    num_unknowns = size(A,2)
    rankAb = rank([A, b])
    dof = num_unknowns - rankAb
end

if ~existence
    error('Error! the forces in getFrictionlessRxn_3d are inconsistent, no solution exists, stopping.');
end

% otherwise, get the psuedoinverse solution.
% If solutions exist, we know this is one of them, which minimizes the
% 2-norm of R.
R = A \ b;

% and parse out into individual blocks.
rx = R(1 : v);
ry = R(v+1 : 2*v);
rz = R(2*v+1 : 3*v);

if( debugging >= 2 )
    disp('Reaction force solutions from getFrictionlessRxn_3d are:');
    rx
    ry
    rz
end
    

end

























% get_reaction_forces_3d.m
% Copyright Andrew P. Sabelhaus 2018

function [px, py, pz] = get_reaction_forces_3d(coordinates, pinned, m, g, debugging)
% get_reaction_forces calculates the external reaction forces on a
% tensegrity structure under the influence of gravity.
%
% This function makes a variety of assumptions about the tensegrity in
% question and its problem setup. Do not use it unless you've ensured those
% assumptions hold (example - distribution of mass, etc.) See various
% papers by Friesen, Sabelhaus, etc. for more discussion.
%   
% Inputs:
%   coordinates = 3 x n matrix of the coordinates of each node.
%   pinned = n x 1 vector, boolean, with 1 = reaction forces exist at that
%       node (i.e., it is pinned) or with = 0 free floating.
%   m = n x 1 vector of the mass of each individual node
%   g = gravitational constant.
%   debugging = verbosity. 0 = none, 1 = messages, 2 = info about matrices
%       and vectors being calculated.
%
% Outputs:
%   px = n x 1 vector of reaction forces in the x direction at each node
%   py = same, for y
%   pz = same, for z
%

% NOTE: this force balance is NOT the same directionality as the cable net
% itself. Gravity is in a different direction.
% For example:
%
%   Cable net solves \sum F_z = 0 at node.
%   Free body diagram has positive q <=> cable forces acting in +
%   direction for each coordinate.
%   this can be seen by doing an example, with C^\top diag(Cz).
%   External forces (p) are assumed to act in the + direction for each
%   coordinate.
%   so \sum cables_z + \sum F_ext_z = 0
%   rearranging, F_ext_z = - F_c_z
%   this is what we get from the configuration matrix, as expected.
%   Since here, F_ext_z = -mg,
%   We'd get a positive tension force q. All OK.
%
%   However: for reaction balance, same thing, but the configuration matrix
%   does NOT take care of the signs for us.
%   Rxn_y + F_ext_z = 0
%   Rxn_y = -F_ext_z
%   Rxn_y = mg
%
%   CONCLUDE: everywhere the external force vector p is used below, we
%   change it to -p!!!!
%   Intuition is that we're balancing out forces, so the reaction force
%   scalars must be the same magnitude as gravity, since their vector
%   directions are opposite.

% Steps:
%   1) calculate sum of forces: x,y are none, z is total gravity (-mg)
%   2) calculate total moment contribution due to point masses
%   3) calculate moment contribution due to (unknown) rxn forces
%   4) solve, using least squares, for rxn forces and plug in.

%% setup

% number of nodes
n = size(coordinates, 2);

% number of nodes that have reaction forces. Since 'pinned' is ones and
% zeros, can just count the number of ones.
v = sum(pinned);

% quick checks
if size(m) ~= n
    error('Error, mass vector should have same number of elements as the number of nodes.');
end

if v > n
    error('Error, the vector of pinned nodes is inconsistent with total number of nodes.');
end

% Take out the coordinates into individual vectors:
% We want them as COLUMNS
x = coordinates(1,:)';
y = coordinates(2,:)';
z = coordinates(3,:)';

% total mass is 
m_tot = sum(m);

% A debugging message if desired, and if the problem is statically
% indeterminate. (One case when this happens is more columns than rows)
if debugging >= 1
    % There are 3*v reaction forces, and 6 balances (Fx,y,z, Mx,y,z)
    if 3*v > 6
        disp('Note, your problem may be statically indeterminate, with respect to reaction forces, since the number of reaction forces');
        disp('is greater than 6, the number of force/moment balances. The solutions for p are then not unique, and it is unknown how this effects the final q.');
        disp('This situation is common and does not indicate that solutions do not exist; however, consider removing the nodes');
        disp('with reaction forces (treat them as anchors.)');
    end
end

%% 1)

% Total forces in x, y, z.
% We're solving Ar R = br.
% Since this is for the entire structure, there's only one sum in each
% dimension and one moment in each dimension, so
br = zeros(6,1);

% The first and second dimensions are zero (sum of forces in x,y is zero.)
% The third dimension is sum of z forces,
F_grav = -m_tot * g;

% NOTE: as discussed above, we must switch the signs here, unlike with the
% force density balance.
br(3) = -F_grav;

% The third element plugged in in section 2 below.

% The constraint matrix Ar is of size (num balances, num rxn forces)
Ar = zeros(6, 3*v);

% since there are 6 sums, and three reaction forces per pinned node.

% Now we plug in for the first three rows of Ar.
% It should look like 
% [1, 1, ... 1, 0, 0, ... 0, 0, 0, ... 0;
%  0, 0, ... 0, 1, 1, ... 1, 0, 0, ... 0;
%  0, 0, ... 0, 0, 0, ... 0, 1, 1, ... 1];

% ...e.g., we're letting the decision variables from 1 to v be the reaction
% forces in x, from v+1 to 2v being the reaction forces in y,
% and 2v+1 to 3v being the reaction forces in z.

% So we just want to sum up the three blocks.
% In other work in this library, we've noted a great way to do this using
% the kronecker product.
% The three top rows of Ar:
Ar(1:3, :) = kron(eye(3), ones(v,1)');


%% 2)

% We're going to sum moments about the origin.

% The RHS of Ar R = br
% Moments in 3D can be expressed as a skew-symmetric matrix, which has the form
% below.

% The external forces at each node are
ext_Fx = zeros(n,1);
ext_Fy = zeros(n,1); 
ext_Fz = m*g;
% ****NOTE**** that we've done the same change as F_grav -> br(3) here,
% that is, changed the sign already so it's +mg.
% (recalling that m is already vector'd out)
ext_F = [ext_Fx; ext_Fy; ext_Fz];

% We'll multiply by a matrix B to get the moment from just this force. 
% We want to evaluate something that looks
% like a cross product, in matrix form of
% [ 0,  -z,   y;
%   z,   0,  -x;
%  -y,   x,   0];

% Since we want the sum of moments around the origin, the moment arm is
% just the coordinate of each node.
% Summing all moments means we output a 3x1 vector. So, taking in
% ext_F \in R^{3n}, means that B should be \in R^{3 x 3n}.
% That makes sense. So, instead of individual moments, we  can use the
% coordinates as whole vectors.

% a quick notational convenience
zos = zeros(n,1);

B = [ zos',  -z',     y';
      z',     zos',  -x';
     -y',     x',     zos'];


% Finally, the (scalar) sum of moments is
ext_M = B * ext_F;

% plug in. These should be the last 3 entries into br.
br(4:end) = ext_M;

                    
%% 3)

% For the actual reaction forces (moment arms),
% first isolate the nodes at which rxn forces exist.

% The construction
% ~any(pinned,2)
% returns a vector of zeros and ones, where the ones are the coordiantes to
% remove.
% So, remove nodes from the x and y vectors.
x_pinned = x;
y_pinned = y;
z_pinned = z;
x_pinned(~any(pinned,2)) = [];
y_pinned(~any(pinned,2)) = [];
z_pinned(~any(pinned,2)) = [];

% The construction is the same as the point masses,
% where here the Fx and Ry are part of the vector R, rxn forces,
% so the row in the constraint matrix is just B but for only the pinned
% nodes.
% This part of Ar is 3 rows x 3v columns, which matches the above for the
% point masses but with v instead of n.

% another quick notational convenience for a vector of zeros for the
% 'pinned' formulation
zosp = zeros(v,1);

B_p = [ zosp',     -z_pinned',   y_pinned';
        z_pinned',  zosp',      -x_pinned';
       -y_pinned',  x_pinned',   zosp'];

% and plug in:
Ar(4:end, :) = B_p;


%% 4) 

if debugging >= 2
    Ar
    br
    rank(Ar)
end

% FINALLY, solve:
R = Ar \ br;

% note, R is \in 3v x 1.

% TO-DO: could actually formulate this as a minimization problem since it's
% statically indeterminate (e.g. nonzero null space, there are *many* more
% rows than columns since 3v >= 6 for almost all applications because at
% least two nodes would be pinned), so could minimize the 2-norm of force.
% Currently unknown if this has any effect on the inverse kinematics
% solution.


%% 5) 

% Now, R is a 3v x 1 vector of px1 ... pxv, py1 ... pyv, pz1 ... pzv
% We need to redistribute each component back into the appropriate index as
% per the 'pinned' vector.
% Let's loop through (again)
% TO-DO: performance improvements, remove looping.
% keep a counter for the node to insert
next_node = 1;

% solutions are
px = zeros(n, 1);
py = zeros(n, 1);
pz = zeros(n, 1);

for i=1:n
    if pinned(i)
        % insert the next set of forces from R into the solutions vector
        % The R vector has the x, y, z, forces stacked in blocks of v.
        px(i) = R(next_node);
        py(i) = R(next_node + v);
        pz(i) = R(next_node + 2*v);
        % increment counter
        next_node = next_node + 1;
    end
end

% a quick check: next_node should now be larger than v.
if next_node <= v
    error('Error in placing node reaction forces into vectors, bug.');
end

% some checking
if debugging >= 2
    R
    px
    py
end
    
end










% getRxn2d.m
% Copyright Andrew P. Sabelhaus 2018

function [px, py] = getRxn2d(x, y, pinned, m, g, debugging)
% getRxn2d calculates the external reaction forces on a
% tensegrity structure under the influence of gravit, in two dimensions.
%
%   ****** NOTE: this function is unused and may still contain bugs. Do not
%           use unless you have vetted it yourself.
%
%   ****** The style in this file has NOT been standardized with the rest
%           of the library. Please excuse the temporary awkwardness.
%
% This function makes a variety of assumptions about the tensegrity in
% question and its problem setup. Do not use it unless you've ensured those
% assumptions hold (example - distribution of mass, etc.) See various
% papers by Friesen, Sabelhaus, etc. for more discussion.
%   
% Inputs:
%   x = component-wise vector of x coordinates for each node (\in R^n)
%   y = component-wise vector of y coordinates for each node (\in R^n)
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
%

% NOTE: this force balance is NOT the same directionality as the cable net
% itself. Gravity is in a different direction.
% For example:
%
%   Cable net solves \sum F_y = 0 at node.
%   Free body diagram has positive q <=> cable forces acting in +
%   direction for each coordinate.
%   this can be seen by doing an example, with C^\top diag(Cy).
%   External forces (p) are assumed to act in the + direction for each
%   coordinate.
%   so \sum cables_y + \sum F_ext_y = 0
%   rearranging, F_ext_y = - F_c_y
%   this is what we get from the configuration matrix, as expected.
%   Since here, F_ext_y = -mg,
%   We'd get a positive tension force q. All OK.
%
%   However: for reaction balance, same thing, but the configuration matrix
%   does NOT take care of the signs for us.
%   Rxn_y + F_ext_y = 0
%   Rxn_y = -F_ext_y
%   Rxn_y = mg
%
%   CONCLUDE: everywhere the external force vector p is used below, we
%   change it to -p!!!!
%   Intuition is that we're balancing out forces, so the reaction force
%   scalars must be the same magnitude as gravity, since their vector
%   directions are opposite.

% Steps:
%   1) calculate sum of forces: x is none, y is total gravity (-mg)
%   2) calculate total moment contribution due to point masses
%   3) calculate moment contribution due to (unknown) rxn forces
%   4) solve, using least squares, for rxn forces and plug in.

%% setup

% number of nodes
n = size(x, 1);

% number of nodes that have reaction forces.
v = nnz(pinned);

% quick checks
if size(m) ~= n
    error('Error, mass vector should have same number of elements as the number of nodes.');
end

% total mass is 
m_tot = sum(m);

if debugging >= 2
    v
    coords_v
end

% A debugging message if desired, and if the problem is statically
% indeterminate. (One case when this happens is more columns than rows)
if debugging >= 1
    % There are 2*v reaction forces, and 3 balances (Fx, Fy, M)
    if 2*v > 3
        disp('Note, your problem may be statically indeterminate, with respect to reaction forces, since the number of reaction forces');
        disp('is greater than 3, the number of force/moment balances. The solutions for p are then not unique, and it is unknown how this effects the final q.');
        disp('This situation is common and does not indicate that solutions do not exist; however, consider removing the nodes');
        disp('with reaction forces (treat them as anchors.)');
    end
end

%% 1)

% Total forces in x and y.
% We're solving Ar R = br.
% Since this is for the entire structure, there's only one sum in each
% dimension and one moment, so
br = zeros(3,1);

% The first dimension is zero (sum of forces in x is zero.)
% The second dimension is sum of y forces,
F_grav = -m_tot * g;

% NOTE: as discussed above, we must switch the signs here, unlike with the
% force density balance.
br(2) = -F_grav;

% The third element plugged in in section 2 below.

% The constraint matrix Ar is of size
Ar = zeros(3, 2*v);

% since there are 3 sums, and two reaction forces per pinned node.

% Now we plug in for the first two rows of Ar.
% It should look like 
% [1, 1, .... 1, 0, 0, ... 0;
%  0, 0, .... 0, 1, 1, .... 1];

% ...e.g., we're letting the decision variables from 1 to v be the reaction
% forces in x, and from v+1 to 2v being the reaction forces in y.
% (note, this changed from original formulation in the 2d horizontal spine
% example, where we grouped according to node instead of according to
% dimension.)

% So we just want to sum up the two blocks.
% In other work in this library, we've noted a great way to do this using
% the kronecker product.
% The two top rows of Ar:
Ar(1:2, :) = kron(eye(2), ones(v,1)');

%% 2)

% The RHS of Ar R = br
% Moments in 2D can be expressed as a determinant, which has a matrix form
% below. This comes out to -yFx + xFy for all nodes i.

% The external forces at each node are
ext_Fx = zeros(n,1);
ext_Fy = m*g;
% ****NOTE**** that we've done the same change as F_grav -> br(2) here,
% that is, changed the sign already so it's +mg.
% (recalling that m is already vector'd out)
ext_F = [ext_Fx; ext_Fy];
% We'll multiply by a matrix
B = [-y', x'];

% Finally, the (scalar) sum of moments is
ext_M = B * ext_F;

% plug in.
br(3) = ext_M;

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
x_pinned(~any(pinned,2)) = [];
y_pinned(~any(pinned,2)) = [];

% The construction is the same as the point masses, -yFx + xFy,
% where here the Fx and Fy are part of the vector R, rxn forces,
% so the row in the constraint matrix is just B but for only the pinned
% nodes.
B_p = [-y_pinned', x_pinned'];

% and plug in:
Ar(3, :) = B_p;

% this is consistent with the done-by-hand example that originally existed
% in horizontal_spine_invkin_2d.
                    
%% 4)

if debugging >= 2
    Ar
    br
    rank(Ar)
end

% FINALLY, solve:
R = Ar \ br;

% note, R is \in 2v x 1.

% TO-DO: could actually formulate this as a minimization problem since it's
% statically indeterminate (e.g. nonzero null space, there are *many* more
% rows than columns since 2v >= 3 for almost all applications because at
% least two nodes would be pinned), so could minimize the 2-norm of force.
% Currently unknown if this has any effect on the inverse kinematics
% solution.


%% 5) 

% Now, R is a 2v x 1 vector of Fx1 ... Fxv, Fy1 ... Fyv]
% We need to redistribute each component back into the appropriate index as
% per the 'pinned' vector.
% Let's loop through.
% TO-DO: performance improvements, remove looping.
% keep a counter for the node to insert
next_node = 1;

% solutions are
px = zeros(n, 1);
py = zeros(n, 1);

for i=1:n
    if pinned(i)
        % insert the next set of forces from R into the solutions vector
        % The R vector has the x, y forces stacked in blocks of v.
        px(i) = R(next_node);
        % The next_node-th position in the second block of R is at index of
        % v plus next_node 
        % (think: 8 nodes = n, 5 nodes pinned = v,
        % insertion of third pinned node (next_node = 3) which is at 4th
        % node overall (i = 4) will be at indices in R of 3 and 8
        py(i) = R(next_node + v);
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










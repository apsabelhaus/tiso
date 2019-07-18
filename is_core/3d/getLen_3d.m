% Andrew P. Sabelhaus 2019

function len = getLen_3d(x, y, z, s, C)
%% getLen_3d
%   getLen_3d calculates the lengths of the cables within a structure, for
%   a three-dimensional problem.
%
%   Inputs:
%       x, y, z = nodal coordinates for the whole structure, \in \mathbb{R}^n
%
%       s = number of cables. Used to create the H matrix that removes the
%       bars from these calculations.
%
%       C = connectivity matrix (also called incidence matrix), \in
%       \mathbb{R}^{(s+r) x n}. See paper for detailed description.
%
%   Outputs:
%
%       len = lengths of each cable. \in \mathbb{R}_+^n.


% validation:
% TO-DO

% The number of total members is the rows of C
m = size(C, 1);
% the number of bars and number of cables have to add up to total members
r = m - s;

% then the H matrix for bar removals can be generated
H = get_H(s, r);

% the component-wise signed lengths of each cable are
dx = H' * C * x;
dy = H' * C * y;
dz = H' * C * z;

% so the lengths of each cable are the euclidean norm of each 2-vector.
% re-organize:
D = [dx, dy, dz];

% the scalar lengths are then the 2-norm (euclidean) for each column, which
% is
len = vecnorm(D, 2, 2);

end


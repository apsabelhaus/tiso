% Andrew P. Sabelhaus 2019

function len = getLen_2d(x, y, H, C)
%% getLen_2d
%   getLen_2d calculates the lengths of the cables within a structure, for
%   a two-dimensional problem.
%
%   Inputs:
%       x, y = nodal coordinates for the whole structure, \in \mathbb{R}^n
%
%       H = bar removal matrix (see get_H), \in \mathbb{R}^{(s+r) x s}.
%
%       C = connectivity matrix (also called incidence matrix), \in
%       \mathbb{R}^{(s+r) x n}. See paper for detailed description.
%
%   Outputs:
%
%       len = lengths of each cable. \in \mathbb{R}_+^n.


% validation:

% The number of nodes, according to the C matrix, is its number of columns:
n = size(C,2);
% Verify that this is also the length of all other nodal vectors.
if n ~= size(x, 1)
    error('Error: the C matrix and the x vector (node positions in x) have a different number of nodes. Cannot continue.');
elseif n ~= size(y, 1)
    error('Error: the C matrix and the y vector (node positions in y) have a different number of nodes. Cannot continue.');
elseif size(x, 2) ~= 1
    error('Error: x is not a column vector. Cannot continue.');
elseif size(y, 2) ~= 1
    error('Error: z is not a column vector. Cannot continue.');
end

% the component-wise signed lengths of each cable are
dx = H' * C * x;
dy = H' * C * y;

% so the lengths of each cable are the euclidean norm of each 2-vector.
% re-organize:
D = [dx, dy];

% the scalar lengths are then the 2-norm (euclidean) for each column, which
% is
len = vecnorm(D, 2, 2);

end


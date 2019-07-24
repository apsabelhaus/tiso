% get_Kb.m
% Copyright Andrew P. Sabelhaus, 2019

function [Kb, K] = get_Kb(eta, d)
%% get_Kb.m
%
%   Generates the K and Kb matrices for a multiple-eta tensegrity
%   structure, i.e., a structure with internal bending moments BUT that has
%   a different number of nodes per body.
%
%   Inputs:
%
%       eta \in \mathbb{Z}_+^(b) = a vector, in each entry is the number of
%       nodes for that body.
%
%       d = {2 or 3}, dimensionality of the problem, for creating K from
%       Kb.
%
%   Outputs:
%
%       Kb = the collapsation matrix for all bodies, in one dimension. This
%       is a diagonal patterning of ones vectors of the size of each 
%       body's eta.
%
%       K = pattern og Kb for each dimension.
%
%   Example:
%
%       eta = [5; 3], d = 3:
%
%       Kb = [1 1 1 1 1, zeros;
%             zeros,     1 1 1];
%
%       K = [Kb, zeros, zeros;
%            zeros, Kb, zeros;
%            zeros, zeros, Kb];
%
%       Note, this reduces to the "same number of nodes in all bodies" case
%       when setting eta = ones(b,1) * (number of nodes in a body.)

% validate: eta should be a column
if ~( size(eta,2) == 1 )
    error('Error: eta must be a column vector.');
end

% the number of bodies is
b = size(eta,1);

% total number of nodes is the sum of all nodes in each body,
n = sum(eta);

% preallocate Kb then, has same rows as number of bodies and one column per
% node (since each node belongs to one body.)
Kb = zeros(b, n);

% Insert the ones vectors into Kb.
% let's keep a running count of how many nodes have been inserted, makes it
% easier to set the columns for the next body.
nn = 0;

for i = 1:b
    % the number of nodes for this body is
    eta_i = eta(i);
    % the indices of the columns here are the next node up till now, then
    % up to adding in the number of nodes for this body
    nn_next = nn + eta_i;
    Kb( i, (nn+1 : nn_next) ) = ones(eta_i,1)';
    % increment count
    nn = nn_next;
end

% pattern out for each dimension.
K = kron(eye(d), Kb);

end













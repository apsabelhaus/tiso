% getObj_2norm.m
% Copyright Andrew P. Sabelhaus, 2019

function weightMat = getObj_2norm(s)
%% getObj_2norm.m
%
%   getObj_2norm calculates the weighting matrix for a quadratic objective
%   function in an optimization problem, with the following form:
%
%   min q_s^\top q_s
%
%   ...so, as x^\top H x, simply H = eye(s).

% validate: s should be a scalar
if ~isscalar(s)
    error('Error: getObj_2norm expected a scalar number of cables s, but s is not scalar.');
end

weightMat = eye(s);

end


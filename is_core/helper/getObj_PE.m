% getObj_PE.m
% Copyright Andrew P. Sabelhaus, 2019

function weightMat = getObj_PE(s, len, kappa)
%% getObj_PE.m
%
%   getObj_PE calculates the weighting matrix for a quadratic objective
%   function in an optimization problem, so that the objective is the total
%   potential energy in the cables. Assumes linear elastic springs.
%   This is of the form:
%
%   min q_s^\top Kappa^\-1 L^s q_s = sum_i kappa_i*(\ell_i - u_i)^2
%
%   Inputs:
%
%       s = number of cables. Used for validating the other two inputs and
%       for consistency with other functions - technically not needed...
%
%       len = vector of lengths of each cable. \in \mathbb{R}_+^s
%
%       kappa = spring constant(s) for the cables. Function is polymorphic
%       with either:
%           \in \mathbb{R}_+ : one spring constant for all cables
%           \in \mathbb{R}_+^s : vector of different spring constants
%                                per-cable
%
%   Outputs:
%
%       weightMat = the 'H' matrix for q_s^\top H q_s.

% validate:
if ~isscalar(s)
    error('Error: getObj_PE expected a scalar number of cables s, but s is not scalar.');
elseif size(len, 1) ~= s
    error('Error: the lengths vector needs to contain all cables, size s.');
elseif size(len, 2) ~= 1
    error('Error: len is not a column vector. Cannot continue.');
elseif any(len < 0)
    error('Error: negative length passed in to getObj_PE.');
end

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

% Expand kappa if a scalar was passed in.
if isscalar(kappa)
    kappa = kappa * ones(s,1);
end

% As per the math in our various papers, total potential energy in cables
% is

kappaInv = 1./kappa;
KappaInv = diag(kappaInv);

lenSqr = len.^2;
LenSqr = diag(lenSqr);

weightMat = KappaInv * LenSqr;

end








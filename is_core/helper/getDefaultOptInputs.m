% Andrew P. Sabelhaus 2019

function oi = getDefaultOptInputs(n, s)
%% getDefaultOptInputs
%
%   returns a struct of the default arguments for the optional struct
%   passed in to an ISO problem.
%
%   Inputs:
%       n = number of nodes
%       s = number of cables
%
%   Outputs:
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

oi = struct();
oi.w = ones(n,1);
oi.R = getObj_2norm(s);
oi.qMin = 0;
oi.debugging = 1;

end


% trajStraightMultiBody_3d.m
% Andrew Sabelhaus 2019

function [xiAll] = trajStraightMultiBody_3d(translation, rotation, b)
%% trajStraightMultiBody_3d
%
%   This function returns a trajectory for a 3d, multi-body tensegrity that 
%   is 'straight', cantilevered out.
%   It includes full position state information, and is a kinematic
%   trajectory (NOT dynamic.)
%
% Inputs:
%
%   translation0 = translation to pattern for CoM of each body. For the
%       horizontal spine, a good choice is barEndpoint * (2/3).
%       Should be \in R^3 for three dimensions.
%
%   rotation0 = rotation angles for all rigid bodies. 
%       Should be \in subset of R^3 for three dimension.
%
%   rotAxisPt = the point about which to apply the rotation for the
%       moving bodies. This is added here to create the possibility for
%       extrinsic rotations around somewhere besides the origin. Column
%       vector, [x; y], \in R^2.
%
%   b = number of bodies to include in the trajectory. Must be => 2.
%
% Outputs:
%   xiAll = system states (kinematic! no velocities!) for each body.
%       That's 6 per body, b bodies, so \xi \in R^{6b}.

% Note that we actually need to have b+1 bodies here.
% TO-DO: make this more extensible, not assuming there's exactly one body
% anchored...
xiAll = zeros(6*(b+1), 1);
% the first vertebra is rotated too
xiAll(4:6, 1) = rotation;

% For each body,
for k = 1:b
    % This body's state is
    xi_k = zeros(6, 1);
    % translate it out by b
    xi_k(1:3) = k * translation;
    % the rotation is the same for each body
    xi_k(4:6) = rotation;
    % plug in to the total state vector
    xiAll(6*k+1 : 6*k+6, 1) = xi_k;
end
    
end

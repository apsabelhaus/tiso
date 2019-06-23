% trajBend2b_2d.m
% Andrew Sabelhaus 2019

% This function returns a trajectory for a 2d, two-body (one moving body, 
% one fixed) tensegrity that bends around the Y+ axis, in either direction.
% It includes full position state information, and is a kinematic
% trajectory (NOT dynamic.)

function [xiAll] = trajBend2b_2d(translation0, rotation0, minSweep, maxSweep, rotAxisPt, numPts)
% Inputs:
%   translation0 = translation of CoM of moving body, at sweep = 0. For the
%       horizontal 2d spine, a good choice is barEndpoint * (2/3).
%       Should be \in R^2 for two dimensions.
%   rotation0 = rotation angle for *both* rigid bodies at sweep = 0. The
%       motivation here is that the local frame can be specified
%       arbitrarily for each body, then orientation can be corrected in
%       this script. (It's just plugged into the return value.)
%   min,maxSweep = the min/max sweep angle for the bend. Trajectory starts at
%       just 'translation' for sweep = minSweep, then rotates around the origin
%       (NOT the local frame) up to maxSweep.
%   rotAxisPt = the point about which to apply the rotation for the
%       moving vertebra. This is added here to create the possibility for
%       extrinsic rotations around somewhere besides the origin. Column
%       vector, [x; y], \in R^2.
%   numPts = the number of timesteps/waypoints in this trajectory.
%
% Outputs:
%   xiAll = system states (kinematic! no velocities!) for each body.
%       That's 3 per body, 2 bodies, so \xi \in R^{6 x numPts}.

% The "first" body is assumed not to move, has its center (or wherever the
% local frame is references from) at the origin, with a rotation given from
% the input.
xiAll = zeros(6, numPts);
% First body rotation. xi is [x, y, theta]
xiAll(3,:) = rotation0;

% We'll span out all the angles to insert.
% beta is sweep angle.
%beta_0 = 0;
beta0 = minSweep;
beta = linspace(beta0, maxSweep, numPts)';

% In order to keep the moving body parallel to the sweept-out centerling
% (e.g. rotating around the origin), the sweep angle and initial angle add.
% abuse some matlab functionality which patterns out scalars to vectors
% automatically:
xiAll(6,:) = rotation0 + beta;

% Finally, do all the positions of the center (or, origin of local frame)
% for the moving body.
for i=1:numPts
    % At the i-th angle, CoM position is translation rotated around the 
    % given axis by this beta. One way to do this is:
    % 1) move frame "out" by axis amount
    % 2) do extrinsic rotation
    % 3) move frame "back" by axis amount
    
    % 1)
    freeCenter = -rotAxisPt + translation0;
    
    % 2)
    xiAll(4:5,i) = [cos(beta(i)),    -sin(beta(i));
                     sin(beta(i)),     cos(beta(i))] * freeCenter;
                 
    % 3)
    xiAll(4:5,i) = xiAll(4:5,i) + rotAxisPt;

    % aren't rotation matrices nice?
end
    
end

% trajBendMultiBody_2d.m
% Andrew Sabelhaus 2019

function [xiAll] = trajBendMultiBody_2d(translation0, rotation0, minSweep, maxSweep, rotAxisPt, numPts, b)
%% trajBendMultiBody_2d
%
%   This function returns a trajectory for a 2d, multi-body tensegrity that 
%   bends around the Y+ axis, in either direction.
%   It includes full position state information, and is a kinematic
%   trajectory (NOT dynamic.)
%
% Inputs:
%
%   translation0 = translation to pattern for CoM of each body, at sweep = 0. For the
%       horizontal 2d spine, a good choice is barEndpoint * (2/3).
%       Should be \in R^2 for two dimensions.
%
%   rotation0 = rotation angle for all rigid bodies at sweep = 0. The
%       motivation here is that the local frame can be specified
%       arbitrarily for each body, then orientation can be corrected in
%       this script. (It's just plugged into the return value.)
%
%   min,maxSweep = the min/max sweep angle for the bend. Trajectory starts at
%       just 'translation' for sweep = minSweep, then rotates around the origin
%       (NOT the local frame) up to maxSweep.
%       These are for the FURTHEST OUT body: all others are 
%
%   rotAxisPt = the point about which to apply the rotation for the
%       moving bodies. This is added here to create the possibility for
%       extrinsic rotations around somewhere besides the origin. Column
%       vector, [x; y], \in R^2.
%
%   numPts = the number of timesteps/waypoints in this trajectory.
%
%   b = number of bodies to include in the trajectory. Must be => 2.
%       THIS INCLUDES THE NON-MOVING BODY. Example, 2 moving bodies, b=3.
%
% Outputs:
%   xiAll = system states (kinematic! no velocities!) for each body.
%       That's 3 per body, b bodies, so \xi \in R^{3b x numPts}.

% The "first" body is assumed not to move, has its center (or wherever the
% local frame is references from) at the origin, with a rotation given from
% the input.
xiAll = zeros(3*b, numPts);
% First body rotation. xi is [x, y, theta]
% xiAll(3,:) = rotation0;

% We'll span out all the angles to insert.
% beta is sweep angle.

% % Ratio of change between each body, as a way to scale the sweep.
% % For example, in the T-CST paper on the spine we did three vertebrae
% % with max sweeps pi/16, pi/12, pi/8, which is roughly 2/3 to 3/4 per
% % body...
% % sweepRatio = 2/3;
% % sweepRatio = 1/2;
% sweepRatio = 3/4;
% 
% % Initialize the sweep angles
% % timesteps are columns
% beta = zeros(b, numPts);
% % for each body, counting backwards from the largest change...
% nextMax = maxSweep;
% nextMin = minSweep;
% % decrementing, but not doing the non-moving body.
% for i=b:-1:2
%     % scaling the min and max sweep for this body.
%     % For example,
%     betaMin_i = nextMin;
%     betaMax_i = nextMax;
%     % insert into trajectory
%     % BUT also add the initial rotation to correct for initial pose.
%     % In order to keep the moving body parallel to the swept-out centerline
%     % (e.g. rotating around the origin), the sweep angle and initial angle add.
% %     beta(i,:) = linspace(betaMin_i, betaMax_i, numPts) + rotation0;
%     beta(i,:) = linspace(betaMin_i, betaMax_i, numPts);
%     % and into the full state trajectory, where the third coordinate is
%     % rotation. Example, 4th body at 3b+2.
%     xiAll(3*i, :) = beta(i, :);
%     % increment for subsequent bodies.
%     nextMax = nextMax * sweepRatio;
%     nextMin = nextMin * sweepRatio;
% end

% ...another way to do it is to just have the sweep evenly divided between
% the bodies. 

% Initialize the sweep angles
% timesteps are columns
beta = zeros(b, numPts);
% for each body, counting backwards from the largest change...
% nextMax = maxSweep;
% nextMin = minSweep;
% decrementing, but not doing the non-moving body.
for i=2:b
    % scaling the min and max sweep for this body.
    % For example, (accounting for the currently present off by one error
    betaMin_i = minSweep*((i-1)/(b-1));
    betaMax_i = maxSweep*((i-1)/(b-1));
    % insert into trajectory
    % BUT also add the initial rotation to correct for initial pose.
    % In order to keep the moving body parallel to the swept-out centerline
    % (e.g. rotating around the origin), the sweep angle and initial angle add.
%     beta(i,:) = linspace(betaMin_i, betaMax_i, numPts) + rotation0;
    beta(i,:) = linspace(betaMin_i, betaMax_i, numPts);
    % and into the full state trajectory, where the third coordinate is
    % rotation. Example, 4th body at 3b+2.
    xiAll(3*i, :) = beta(i, :);
    % increment for subsequent bodies.
%     nextMax = nextMax * sweepRatio;
%     nextMin = nextMin * sweepRatio;
end

% Finally, do all the positions of the center (or, origin of local frame)
% for the moving body.
for i=1:numPts
    % At the i-th angle, CoM position is translation rotated around the 
    % given axis by this beta. One way to do this is:
    % 1) move frame "out" by axis amount
    % 2) do extrinsic rotation
    % 3) move frame "back" by axis amount
    
    % per body, only doing the moving ones,
    for k=2:b
        % 1)
        % the translation of this body in its initial pose is
        % scaled by how many bodies "out" from the rotation point it will
        % be
        freeCenter_k = -rotAxisPt + translation0*(k-1);
        % 2)
        % its rotated position vector at this timestep is
        pos_ki = [cos(beta(k,i)),    -sin(beta(k,i));
                  sin(beta(k,i)),     cos(beta(k,i))] * freeCenter_k;
        %
        % 3) translating it back
%         xiAll(4:5,i) = xiAll(4:5,i) + rotAxisPt;
        pos_ki_adj = pos_ki + rotAxisPt;
        % 4) insert into state.
        x_ki_index = 3*(k-1) + 1;
        y_ki_index = 3*(k-1) + 2;
        xiAll(x_ki_index:y_ki_index, i) = pos_ki_adj;
    end

    % aren't rotation matrices nice?
end

% Finally, we actually need to add back in the body rotation for each,
% since it was kept separate above.
for k=1:b
    % add rotation0 to this body's angle state.
    xiAll(3*k, :) = xiAll(3*k, :) + rotation0;
end
    
end

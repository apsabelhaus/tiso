% trajectory_XZG_bend_2.m

% This function returns a trajectory for a 2d, single-body tensegrity that
% bends around the Y+ axis, according to the inverse kinematics script, in either direction.
% It includes full position state information, and is a kinematic
% trajectory (NOT dynamic.)

function [xi_all] = trajectory_XZT_bend_2d(translation_0, rotation_0, min_sweep, max_sweep, rot_axis_pt, num_points)
% Inputs:
%   translation_0 = translation of CoM of moving body, at sweep = 0. For the
%       horizontal 2d spine, a good choice is bar_endpoint * (2/3).
%       Should be \in R^2 for two dimensions.
%   rotation_0 = rotation angle for *both* rigid bodies at sweep = 0. The
%       motivation here is that the local frame can be specified
%       arbitrarily for each body, then orientation can be corrected in
%       this script. (It's just plugged into the return value.)
%   max_sweep = the maximum sweep angle for the bend. Trajectory starts at
%       just 'translation' for sweep = 0, then rotates around the origin
%       (NOT the local frame) up to max_sweep.
%   rot_axis_pt = the point about which to apply the rotation for the
%       moving vertebra. This is added here to create the possibility for
%       extrinsic rotations around somewhere besides the origin. Column
%       vector, [x; y], \in R^2.
%   num_points = the number of timesteps/waypoints in this trajectory.
% Outputs:
%   xi_all = system states (kinematic! no velocities!) for each body.
%       That's 3 per body, 2 bodies, so \xi \in R^{6 x num_points}.

% The "first" body is assumed not to move, has its center (or wherever the
% local frame is references from) at the origin, with a rotation given from
% the input.
xi_all = zeros(6, num_points);
% First body rotation. xi is [x, y, theta]
xi_all(3,:) = rotation_0;

% We'll span out all the angles to insert.
% beta is sweep angle.
%beta_0 = 0;
beta_0 = min_sweep;
beta = linspace(beta_0, max_sweep, num_points)';

% In order to keep the moving body parallel to the sweept-out centerling
% (e.g. rotating around the origin), the sweep angle and initial angle add.
% abuse some matlab functionality which patterns out scalars to vectors
% automatically:
xi_all(6,:) = rotation_0 + beta;

% Finally, do all the positions of the center (or, origin of local frame)
% for the moving body.
for i=1:num_points
    % At the i-th angle, CoM position is translation rotated around the 
    % given axis by this beta. One way to do this is:
    % 1) move frame "out" by axis amount
    % 2) do extrinsic rotation
    % 3) move frame "back" by axis amount
    
    % 1)
    free_center = -rot_axis_pt + translation_0;
    % 2)
%     xi_all(4:5,i) = [cos(beta(i)),    -sin(beta(i));
%                      sin(beta(i)),     cos(beta(i))] * translation_0;
    xi_all(4:5,i) = [cos(beta(i)),    -sin(beta(i));
                     sin(beta(i)),     cos(beta(i))] * free_center;
    % 3)
    xi_all(4:5,i) = xi_all(4:5,i) + rot_axis_pt;

    % aren't rotation matrices nice?
end
    
end

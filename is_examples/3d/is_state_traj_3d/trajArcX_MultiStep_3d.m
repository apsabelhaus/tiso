% trajArcX_Multistep_3d.m
% Andrew Sabelhaus 2019

function [xiAll] = trajArcX_MultiStep_3d(x0, xf, h_t0, h_tT, b, numPts)
%% trajArcX_MultiStep_3d
%
%   This function returns a trajectory for a 3d, multi-body tensegrity that 
%   is a parabola in the X-Z plane, Y=0 everywhere. The first body
%   is at [t0, 0, 0] and last body is at [tf, 0, 0], with remaining b-2
%   bodies evenly spaced along the x-axis. Rotations of bodies align with arc
%   tangent points.
%
%   This "multistep" version of the function returns a time series of
%   poses by calling trajArc_X iteratively.
%
% Inputs:
%
%   x0 = x-coordinate of first body
%
%   xf = x-coordinate of last body
%
%   h = apex of the arc (vertex in Z for the parabola)
%       h_t0 = arc apex for first timestep
%       h_tT = arc apex for last timestep
%       * Trajectory is linearly spaced between t0 and tT.
%
%   b = number of bodies to include in the trajectory. Must be => 2.
%
%   numPts = size of the trajectory. Number of timesteps.
%
% Outputs:
%   xiAll = system states (kinematic! no velocities!) for each body.
%       That's 6 per body, b bodies, so \xi \in R^{6b}.

% First, linear spacing over the apex (height of the arc.
arc = linspace(h_t0, h_tT, numPts);

% Iterate over trajArcX_3d and concatenate the results.
xiAll = [];
for i=1:size(arc, 2)
    xiAll(:,end+1) = trajArcX_3d(x0, xf, arc(i), b);
end

end

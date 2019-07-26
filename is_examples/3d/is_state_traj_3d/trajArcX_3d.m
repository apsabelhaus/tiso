% trajArcX_3d.m
% Andrew Sabelhaus 2019

function [xiAll] = trajArcX_3d(x0, xf, h, b)
%% trajArcX_3d
%
%   This function returns a trajectory for a 3d, multi-body tensegrity that 
%   is a parabola in the X-Z plane, Y=0 everywhere. The first body
%   is at [t0, 0, 0] and last body is at [tf, 0, 0], with remaining b-2
%   bodies evenly spaced along the x-axis. Rotations of bodies align with arc
%   tangent points.
%
% Inputs:
%
%   x0 = x-coordinate of first body
%
%   xf = x-coordinate of last body
%
%   h = apex of the arc (vertex in Z for the parabola)
%
%   rotation0 = rotation angles for all rigid bodies. 
%       Should be \in subset of R^3 for three dimension.
%
%   b = number of bodies to include in the trajectory. Must be => 2.
%
% Outputs:
%   xiAll = system states (kinematic! no velocities!) for each body.
%       That's 6 per body, b bodies, so \xi \in R^{6b}.

% There are b bodies, so
xiAll = zeros(6*b, 1);

% Calculate the parameters for the parabola.
% The vertex in the x-direction is halfway between min and max
v = (xf-x0)/2;
% The focus length is a bit more complicated. From Drew's notes,
f = - 1/(4*h) * (x0 - v)^2;
% The equation of the parabola is then 
% z = 1/(4*f) * (x - v)^2 + h

% TEST: just plot the curve.
% num_pts = 100;
% x = linspace(x0, xf, num_pts);
% z = 1/(4*f) * (x-v).^2 + h;
% figure;
% hold on;
% plot(x, z);

% The x-coordinate of each body is evenly spaced
x_coord = linspace(x0, xf, b);
% The z-coordinate can be found via the parabola's curve
z_coord = 1/(4*f) * (x_coord-v).^2 + h;
% The rotations are so the bodies are tangent to the curve, so just the
% derivative of the parabola, rotation is 
rot_y = - 1/(2*f) * (x_coord - v);

% For each body,
for k = 1:b
    % This body's state is
    xi_k = zeros(6, 1);
    % translations in x, z, and rotation in y (assuming YPR angles)
    xi_k(1) = x_coord(k);
    xi_k(3) = z_coord(k);
    xi_k(5) = rot_y(k);
    % plug in to the total state vector
    xiAll(6*(k-1)+1 : 6*(k-1)+6, 1) = xi_k;
end
    
end

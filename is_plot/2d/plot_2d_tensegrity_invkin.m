% plot_2d_tensegrity_invkin.m
% Copyright Andrew P. Sabelhaus 2018

function handles = plot_2d_tensegrity_invkin( C, x, z, s, rad)
%% plot_2d_tensegrity_invkin
%
%   This function plots a 2-dimensional tensegrity structure, with some
%   additional labelling so that it's appropriate for debugging /
%   visualizing an inverse kinematics problem for a tensegrity robot.
%
%   The function creates a new figure window and plots there. The plot is
%   in 3D even though it's a 2D tensegrity, this gives a nice effect and is
%   prettier.
%
%   Technically, MATLAB does "y" as "up", inconsistent with how we usually
%   do robotics, so this script changes "z" to "y". In other words: the
%   coordinate system changes between the calculations and this plot (with
%   respect to MATLAB's assumptions) but this plot should be correct with
%   respect to right-hand coordinate systems and the pictures in our
%   papers.
%
%   Inputs:
%       C = configuration matrix for this tensegrity. This should be for
%       the whole robot, not just one rigid body.
%
%       x, z = locations of each node in the robot.
%
%       s = number of cables in the robot. Since we assume that C is
%       organized by cables-first, the parameter "s" is needed so that the
%       visualization can display an edge as either a cable or bar
%       (depending on where it is in C.)
%
%       rad = radius of the "bars" that we'll plot. This function makes a
%       nice surface illustration for the bar elements, and rad determines
%       how thick the bars are (their radius.)
%
%   Outputs:
%
%       handles = graphics handles to all the objects that are plotted.
%       This is so that the calling function can delete everything if
%       desired. Cell array.
%
%   Depends:
%
%       get_2d_surface_points, a function that calculates the inputs for
%       surf to make a nice bar.
%


%% Setup the problem

% New figure window.
figure;
hold on;

% Pick out the number of nodes (for looping later). This is the number of
% columns in the configuration matrix.
n = size(C, 2);

% similarly, the number of bars is the number of rows minus the number of
% cables.
r = size(C, 1) - s;

% Some constants for the plotting.

% The surfaces are discretized. Let's specify the (amount of discretization
% of the surfaces?) Default is 20?
surf_discretization = 20;
% For the cylinders, we need another discretization for the *length* of the
% cylinder in addition to the *arc length* of the surrounding circles.
surf_length_discretization = 40;

% Specify the color and thickness of cables.
% TO-DO: pass this in?
cable_color = 'r';
cable_thickness = 2;

% Color for the bars:
black = [0,0,0];
bar_color = black;

% The handles array can be a cell array:
handles = {};

% We want the nodes to be a bit bigger than the bars.
% ...but not that much.
node_rad = 1.1 * rad;

% Labeling the nodes needs to happen with some offset from the point
% itself, so that the sphere doesn't overtake the point.
% We'll do the radius of the nodes plus some constant
label_offset = node_rad + 0.02;

% We can also set the color and size of the text.
label_color = 'k';
label_size = 14;

%% Plot the nodes

% A set of matrices for surf-ing a sphere. These will be moved around to
% plot the nodes.
[sphere_x, sphere_y, sphere_z] = sphere(surf_discretization);

% Each sphere will be at rad*(output of sphere) + position offset. 
x_sphere_outer = node_rad * sphere_x;
y_sphere_outer = node_rad * sphere_y;
z_sphere_outer = node_rad * sphere_z;

% Plot spheres at each node.
for i=1:n
    % Translate the sphere positions for surf.
    % NOTE that as described above, we switch y and z since MATLAB doesn't
    % plot using a usual right-handed coordinate system ("y" is "up", here.)
    x_translated = x_sphere_outer + x(i);
    y_translated = y_sphere_outer + z(i);
    z_translated = z_sphere_outer;
    % Plot the surface
    handles{end+1} = surf(gca, x_translated, y_translated, ...
        z_translated, 'LineStyle', 'none', 'edgecolor', bar_color, ...
        'facecolor', bar_color);
    % Put a label for this node.
    % including the offset so it's easier to see.
    handles{end+1} = text(x(i) + label_offset, z(i) + label_offset, 0, ...
        num2str(i), 'Color', label_color, 'FontSize', label_size);
end


%% Plot the cables

% the first s rows of C.
for j=1:s
    % The start and end nodes of this cable are the +1 and -1 entries in
    % this row of C.
    % A neat MATLAB trick here is to compare a vector with 1 or -1, to get
    % a true/false vector, then the 'find' command returns the index of the
    % 'true' element.
    from_index = find( C(j,:) == 1 );
    to_index = find( C(j,:) == -1 );
    % The 'line' command operates on row vectors.
    cable_x = [x(from_index), x(to_index)];
    cable_z = [z(from_index), z(to_index)];
    % ...and also returns a handle that we should be storing.
    handles{end+1} = line(cable_x, cable_z, 'Color', cable_color, ...
        'LineWidth', cable_thickness);
end


%% Plot the bars

% the last r rows of C.
% Example, if there are 4 cables and 6 bars, this should be rows 5 through
% 10 of C. 

for j=(s+1):(s+r)
    % As with the cables, pick out the nodes that this bar will connect:
    % A neat MATLAB trick here is to compare a vector with 1 or -1, to get
    % a true/false vector, then the 'find' command returns the index of the
    % 'true' element.
    from_index = find( C(j,:) == 1 );
    to_index = find( C(j,:) == -1 );
    % The function get_2d_surface_points takes in each point as a column
    % vector \in R^3, so:
    % (and remember we're switching y and z because MATLAB)
    start_pt = [x(from_index), z(from_index), 0];
    end_pt = [x(to_index), z(to_index), 0];
    % Get the inputs for surf-ing this bar
    [x_cyl, y_cyl, z_cyl] = get_2d_surface_points(rad, surf_discretization, ...
                surf_length_discretization, start_pt, end_pt);
    % Finally, plot the bar.
    handles{end+1} = surf(gca, x_cyl, y_cyl, ...
        z_cyl, 'LineStyle', 'none', 'edgecolor', bar_color, ...
        'facecolor', bar_color);
end


%% Some labels
title('Tensegrity Structure used w/InvKin');
xlabel('x position (m)');
ylabel('z position (m)');
% Do a font size change.
set(gca, 'FontSize', 14);

%% Cleanup

% Sometimes this scales funny. Reset the axes.
axis equal;

end















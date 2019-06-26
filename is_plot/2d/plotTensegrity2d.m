% plotTensegrity2d.m
% Copyright Andrew P. Sabelhaus 2018

function handles = plotTensegrity2d( C, x, z, s, rad)
%% plotTensegrity2d
%
%   This function plots a 2-dimensional tensegrity structure, with some
%   additional labelling so that it's appropriate for debugging /
%   visualizing an inverse statics problem for a tensegrity robot.
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
%       getSurfPts2d, a function that calculates the inputs for
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
surfDiscr = 20;
% For the cylinders, we need another discretization for the *length* of the
% cylinder in addition to the *arc length* of the surrounding circles.
surfLengthDiscr = 40;

% Specify the color and thickness of cables.
% TO-DO: pass this in?
cableColor = 'r';
cableThickness = 2;

% Color the nodes differently for illustration.
nodeColor = 'k';

% Color for the bars:
black = [0,0,0];
barColor = black;

% The handles array can be a cell array:
handles = {};

% We want the nodes to be a bit bigger than the bars.
% ...but not that much.
% nodeRad = 1.1 * rad;
nodeRad = 1.2 * rad;

% Labeling the nodes needs to happen with some offset from the point
% itself, so that the sphere doesn't overtake the point.
% We'll do the radius of the nodes plus some constant
labelOffset = nodeRad + 0.02;

% We can also set the color and size of the text.
labelColor = 'k';
labelSize = 14;

%% Plot the nodes

% A set of matrices for surf-ing a sphere. These will be moved around to
% plot the nodes.
[sphereX, sphereY, sphereZ] = sphere(surfDiscr);

% Each sphere will be at rad*(output of sphere) + position offset. 
sphereOuterX = nodeRad * sphereX;
sphereOuterY = nodeRad * sphereY;
sphereOuterZ = nodeRad * sphereZ;

% Plot spheres at each node.
for i=1:n
    % Translate the sphere positions for surf.
    % NOTE that as described above, we switch y and z since MATLAB doesn't
    % plot using a usual right-handed coordinate system ("y" is "up", here.)
    translatedX = sphereOuterX + x(i);
    translatedY = sphereOuterY + z(i);
    translatedZ = sphereOuterZ;
    % Plot the surface
    handles{end+1} = surf(gca, translatedX, translatedY, ...
        translatedZ, 'LineStyle', 'none', 'edgecolor', nodeColor, ...
        'facecolor', nodeColor);
    % ...was barColor for both.
    % Put a label for this node.
    % including the offset so it's easier to see.
%     handles{end+1} = text(x(i) + labelOffset, z(i) + labelOffset, 0, ...
%         num2str(i), 'Color', labelColor, 'FontSize', labelSize);
end


%% Plot the cables

% the first s rows of C.
for j=1:s
    % The start and end nodes of this cable are the +1 and -1 entries in
    % this row of C.
    % A neat MATLAB trick here is to compare a vector with 1 or -1, to get
    % a true/false vector, then the 'find' command returns the index of the
    % 'true' element.
    fromIndex = find( C(j,:) == 1 );
    toIndex = find( C(j,:) == -1 );
    % The 'line' command operates on row vectors.
    cableX = [x(fromIndex), x(toIndex)];
    cableZ = [z(fromIndex), z(toIndex)];
    % ...and also returns a handle that we should be storing.
    handles{end+1} = line(cableX, cableZ, 'Color', cableColor, ...
        'LineWidth', cableThickness);
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
    fromIndex = find( C(j,:) == 1 );
    toIndex = find( C(j,:) == -1 );
    % The function get_2d_surface_points takes in each point as a column
    % vector \in R^3, so:
    % (and remember we're switching y and z because MATLAB)
    startPt = [x(fromIndex), z(fromIndex), 0];
    endPt = [x(toIndex), z(toIndex), 0];
    % Get the inputs for surf-ing this bar
    [cylX, cylY, cylZ] = getSurfPts2d(rad, surfDiscr, ...
                surfLengthDiscr, startPt, endPt);
    % Finally, plot the bar.
    handles{end+1} = surf(gca, cylX, cylY, ...
        cylZ, 'LineStyle', 'none', 'edgecolor', barColor, ...
        'facecolor', barColor);
end


%% Some labels
title('Tensegrity Structure, 2D');
xlabel('X Position (m)');
ylabel('Y Position (m)');
% Do a font size change.
set(gca, 'FontSize', 14);

%% Cleanup

% Sometimes this scales funny. Reset the axes.
axis equal;

end















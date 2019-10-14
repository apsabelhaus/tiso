% plotCableTensions.m
% Copyright Andrew P. Sabelhaus 2019

function hdls = plotCableTensions(fOpt, sigma, maxF, labels)
%% plotCableTensions
%   A helper function that plots a set of cable tensions over time.
%
%   Inputs:
%       fOpt = cable tensions as an s x t matrix, with s cables and t
%       timesteps.
%
%       sigma = how many cables there are per pair of bodies. For example,
%           if there are 4 bodies, connected in a row, probably 3 sets of
%           cables. It's assummed that fOpt is in blocks of size sigma.
%
%       maxF = upper limit on the Y-axis (force).
%
%       labels = a cell array of size sigma, containing strings. One label per each cable in a
%           set. Example, for a spine, these might be "top", "saddle left",
%           etc.
%
%       hdls = a vector of handles to all the plotted data.   
%
%   Outputs: none.


%% Setup the problem

% Adding the gridLegend function here...
addpath( genpath('gridLegend_v1.4') );

% Some constants for the plotting.

% We'll cycle through the markers according to cable within a set,
markers = {'+','o','*','x','v','d','^','s','>','<'};
% and cycle through colors according to pair between bodies.
% colors = {'y', 'm', 'c', 'r', 'g', 'b', 'k'};
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k'};

% line and marker size
lineWidth = 1.5;
markerSize = 6;


% a quick check: this function currently can't handle more than
% size(markers,2) cables in a set
if sigma > size(markers,2)
    error('Error, plotCableTensions cannot currently handle this number of cables in a set between bodies.');
end

% Calculate how many pairs of cables there are
s = size(fOpt,1);
numPairs = s/sigma; % check that this is a whole number...

% set up the window
fontsize = 14;
fig = figure;
hold on;
% Set up the window
set(gca, 'FontSize', fontsize);
% set(fig,'Position',[100,100,500,350]);
set(fig,'Position',[100,100,700,450]);
% set(fig,'PaperPosition',[1,1,5.8,3.5]);
% Annotate the plot
% title('Tensegrity ISO Cable Tensions');
% title('Inverse Statics Opt. Cable Tensions');
title('Inverse Statics Opt. Cable Tensions: Quadruped');
ylabel('Force (N)');
xlabel('Timestep (Pose)');
% legend('Test (Computer Vision)', 'Predicted State', 'Location', 'Best');

% Need an x-axis
T = size(fOpt, 2);
t = 1:T;

% Save the handles to all the data.
hdls = [];

%% Plot per set of cables
for i=1:numPairs
    % the color to use for this pair is
    color_i = colors{i};
    % For each cable in this set, which starts at
    % sigma*(i-1) + 1
    % and ends at
    % sigma*i 
    for k = 1:sigma
        % the marker to use is
        marker_ik = markers{k};
        % and the index to grab from within fOpt is
        cable_ik = sigma*(i-1) + k;
        % so finally, plot it.
        hdls(end+1) = plot(t, fOpt(cable_ik, :), 'Color', color_i, 'Marker', marker_ik, ...
            'MarkerSize', markerSize, 'LineWidth', lineWidth, 'LineStyle', 'none');
    end
end

%% Adjust the plot

% some manual hacks for the plot... for the RA-L paper.
% TO-DO: add arguments or return handles somehow!!

% set(gca, 'FontSize', 16);

ylim([-1, maxF]);
xlim([1, T]);
% TO-DO: ADD LEGEND LOCATION AS AN ARGUMENT!!
[~, objh] = legend(labels, 'Location', 'NW');

% MANUAL HACK FOR 2D Spine
% [~, objh] = legend(labels, 'FontSize', 14, 'Orientation', 'horizontal', 'Location', 'north');

% ANOTHER MANUAL HACK FOR the 3D PLOTS... TO-DO, FIX THIS!!
% numCol3d = 2;
% legend_hdls = gridLegend(hdls(1:sigma), numCol3d, labels, 'TextColor', 'k', 'Location', 'NW');
% child_legend_hdls = get(legend_hdls, 'children');

% Adjust the properties of the legend lines
linesh = findobj(objh, 'type', 'line');
% linesh = findobj(child_legend_hdls, 'type', 'line')
set(linesh, 'Color', 'k', 'LineStyle', 'none');

% hold off;

% % manually adjust the legend font size
% texth = findobj(objh, 'type', 'Text');
% set(texth, 'FontSize', 14);

end















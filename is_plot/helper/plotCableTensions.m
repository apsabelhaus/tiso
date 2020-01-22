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

%%%%%%%%%%%%%%% HARD-CODED FOR RA-L PAPER QUADRUPED
% For the relevant markers:
% For the quadruped:
% alphaMore = 0.7;  % 0.5
% For the mult-vert spine:
alphaMore = 0.6;
% for the ones to de-emphasize:
alphaLess = 0.3; % 0.1

% We'll cycle through the markers according to cable within a set,
% markers = {'+','o','*','x','v','d','^','s','>','<'};

% Used for the initial revision of the RA-L paper, quadruped:
% markers = {'+','o','^','s','.','^','x','>','<','d'};
% For the second revision: are the left, right, and saddle cables even relevant?
% QUADRUPED:
% markers = {'s','o','d','d','.','.','.','.','.','.'};
% 2D MULTIVERT SPINE:
markers = {'s','o','v','^','.','.','.','.','.','.'};

% and cycle through colors according to pair between bodies.
% colors = {'y', 'm', 'c', 'r', 'g', 'b', 'k'};
colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k'};

% line and marker size
lineWidth = 1.5;
% markerSize = 6;
% For the RA-L initial submission:
% markerSize = 10;

% FOR THE RA-L REVISIONS: quadruped
% Now using scatter, the sizes are very different:
markerSize = 140;
% HACK FOR RA-L: adjust the marker size for the irrelevant lines.
% markerSize_smaller = 150; % For scatter: was 100
% For plot, is smaller
% markerSize_smaller = 8;


% a quick check: this function currently can't handle more than
% size(markers,2) cables in a set
if sigma > size(markers,2)
    error('Error, plotCableTensions cannot currently handle this number of cables in a set between bodies.');
end

% Calculate how many pairs of cables there are
s = size(fOpt,1);
numPairs = s/sigma; % check that this is a whole number...

% set up the window
% fontsize = 14;
fontsize = 16;
fig = figure;
hold on;
% Set up the window
set(gca, 'FontSize', fontsize);
% set(fig,'Position',[100,100,500,350]);
set(fig,'Position',[100,100,600,400]);
% set(fig,'Position',[100,100,700,450]);
% set(fig,'PaperPosition',[1,1,5.8,3.5]);
% Annotate the plot
% title('Tensegrity ISO Cable Tensions');
% title('Inverse Statics Opt. Cable Tensions');
% title('Inverse Statics Opt. Cable Tensions: Quadruped');
title('Inverse Statics Opt. Cable Tensions: 2D Spine');
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
        % Two more hacks for the RA-L paper:
        % (1) Don't plot the horizontal left/right cables, they're at qMin
        % and so are not relevant,
        % (2) plot the saddle cables as lines, not markers, so they are
        % easier to see.
        % To do this, here's a variable determining the type of plot to
        % make. 0 = no plot, 1 = marker, 2 = line.
        % FOR THE RA-L QUADRUPED: default will be marker now
        whatToPlot = 1;
        % Hack for RA-L paper: don't plot the horizontal left/right cables,
        % since they're all tensioned at zero.
        % For the RA-L REVISONS QUADRUPED:
%         range_leftright = (3:4);
%         if any( range_leftright == k )
%             % Turn these ones off
%             whatToPlot = 0;
% %             disp('Turning off plotting for k=')
% %             disp(num2str(k))
%         end
        % HACK FOR RA-L PAPER:
        % For the saddle cables, indices 5-8,
        % use smaller markers and more transparency:
        markerSize_adjusted = markerSize;
        markerAlpha_adjusted = alphaMore;
        % FOR RA-L QUADRUPED:
%         testrange = (5:8);
%         % To check if k is in 3:8, we see if any number in 3:8 equals k
%         if any( testrange == k )
%             markerSize_adjusted = markerSize_smaller;
%             markerAlpha_adjusted = alphaLess;
%             % also set these to lines
% %             disp('Doing lines for k=')
% %             disp(num2str(k))
%             whatToPlot = 2;
%         end
        % so finally, plot it.
        if whatToPlot == 1
            thishandle = scatter(t, fOpt(cable_ik, :), markerSize_adjusted, marker_ik, ...
            'MarkerFaceColor', color_i, 'MarkerEdgeColor', color_i);
            % HACK FOR RA-L PAPER: adjust the alpha valua
            thishandle.MarkerFaceAlpha = markerAlpha_adjusted;
            thishandle.MarkerEdgeAlpha = markerAlpha_adjusted;
            hdls(end+1) = thishandle;
        elseif whatToPlot == 2
            thishandle = plot(t, fOpt(cable_ik, :), 'Color', color_i, 'Marker', marker_ik, ...
                'MarkerSize', markerSize_adjusted, 'MarkerFaceColor', color_i, 'LineWidth', lineWidth, 'LineStyle', 'none');
            % FOR RA-L QUADRUPED:
%             thishandle = plot(t, fOpt(cable_ik, :), 'Color', 'k', 'LineWidth', lineWidth, 'LineStyle', '--');
%             thishandle.Color(4) = markerAlpha_adjusted;
            % For the RA-L multivert spine:
            alpha(thishandle, markerAlpha_adjusted);
            hdls(end+1) = thishandle;
        end
        
    end
end

%% Adjust the plot

% some manual hacks for the plot... for the RA-L paper.
% TO-DO: add arguments or return handles somehow!!

% set(gca, 'FontSize', 16);

% Adjust the limits according to size of the markers.
% HACK ON 2019-10-29. To-do: figure out algorithmic way of shifting the
% plot.
% ylim([-1, maxF]);
% xlim([1, T]);
% For the original RA-L paper:
% ylim([-2, maxF-2]);
xlim([0.5, T]);
% For the revision:
% ylim([2, maxF-2]);
% xlim([0.0, T]);
% REVISION, 2D SPINE:
ylim([-3, maxF+2]);



% TO-DO: ADD LEGEND LOCATION AS AN ARGUMENT!!
% Hack for RA-L paper: legend for only the top bottom and saddles
% labels_noleftright = {labels{1:2}, labels{5:8}};
% [~, objh] = legend(labels, 'FontSize', 11, 'Location', 'NW');
% FOR THE RA-L REVISIONS quadruped:
% [~, objh] = legend(labels_noleftright, 'FontSize', 14, 'Location', 'NW');

% MANUAL HACK FOR 2D Spine
[~, objh] = legend(labels, 'FontSize', 14, 'Orientation', 'horizontal', 'Location', 'north');

% ANOTHER MANUAL HACK FOR the 3D PLOTS... TO-DO, FIX THIS!!
% numCol3d = 2;
% legend_hdls = gridLegend(hdls(1:sigma), numCol3d, labels, 'TextColor', 'k', 'Location', 'NW');
% child_legend_hdls = get(legend_hdls, 'children');

% Adjust the properties of the legend lines
linesh = findobj(objh, 'type', 'line');
% For the RA-L Revisions: legends of scatter plots now have 'patch' objects
patchesh = findobj(objh, 'type', 'patch');
% linesh = findobj(child_legend_hdls, 'type', 'line');
% set(linesh, 'Color', 'k', 'LineStyle', 'none');
% FOR THE RA-L QUADRUPED:
% set(linesh, 'Color', 'k');
% FOR THE RA-L MULTI-VERT SPINE:
set(linesh, 'Color', 'k');
% For the RA-L revisions: adjust the face and edge colors of the markers
set(patchesh, 'MarkerFaceColor', 'k', 'LineStyle', 'none');
set(patchesh, 'MarkerEdgeColor', 'k', 'LineStyle', 'none');

% hold off;

% % manually adjust the legend font size
% texth = findobj(objh, 'type', 'Text');
% set(texth, 'FontSize', 14);

end















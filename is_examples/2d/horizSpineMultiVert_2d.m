%% horizSpineMultiVert_2d.m
% Copyright Andrew P. Sabelhaus 2019

% This script used the tiso libraries to solve the inverse statics problem
% for a multi-body (b) tensegrity spine, defined in d=2 dimensions. 
% As the term is used here, 'spine'
% refers to a structure with repeated rigid bodies of the same geometry and

%% set up the workspace
clear all;
close all;
clc;

% add the core libraries
addpath( genpath('../../is_core/2d') );
% and the helper functions
addpath( genpath('../../is_core/helper') );
% same for the plotting
addpath( genpath('../../is_plot/2d') );
addpath( genpath('../../is_plot/helper') );
% the trajectories are now in a subfolder
addpath( genpath('is_state_traj_2d') );

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 1;

% If appropriate, output a starting message.
if debugging >= 1
    disp('Starting horizontal spine 2d rigid body inverse statics example...');
end

% minimum cable force density
%qMin = 0; % units of N/m, depending on m and g
qMin = 0.5;

% For the rest length constraint,
uMin = 0.01;

% Local frame for one rigid body (locations of nodes)

% We'll use a similar frame as Drew's 2018 T-CST paper for a single 2D
% spine. Nodes are column vectors [x; z], but for ease, written as transposed
% here.

%barEndpoint = 0.5; % meters. 50 cm.

% For the 2D spine control test, fall 2018, the dimensions of the vertebra
% are as follows:
% barEndpoint = 4 inches
% that's 4 * 2.54 * 0.01 = 0.1016 meters
barEndpoint = 0.1016;

% local frame of nodes
a = [   0,              0;
        barEndpoint,   -barEndpoint;
        -barEndpoint,  -barEndpoint;
        0,              barEndpoint]';
    

% (e.g., vertebra is in an 8x8 inch box.)

% Mass of a single vertebra
% Mass as measured with a scale on 2018-11-18 is about 500g
mBody = 0.495;
% mBody = 1;
    
if debugging >= 2
    a
end

% number of rigid bodies
% b = 3 including the non-moving one, but since we're removing it later,...
% b = 2;

b = 4;
% When removing the anchor nodes, it's like removing one of the bodies:
% b = 1;

% Configuration matrix. Split into source, sink, and body matrices for each
% vertebra,

% cable source matrix
Cso = [0 1 0 0;
       0 0 1 0;
       0 0 0 1;
       0 0 0 1];

% cable sink matrix
Csi = [0 -1 0 0;
       0 0 -1 0;
       0 -1 0 0;
       0 0 -1 0];

% vertebra rod matrix
Crb = [1 -1 0 0;
       1 0 -1 0;
       1 0 0 -1];
   
% Get the full connectivity matrix.
% NOTE HERE that though we'll be removing the leftmost vertebra later, we
% still need its connectivity matrix blocks, so number of bodies is
% actually b+1.

C = get_C(Cso, Csi, Crb, b+1);

% % Rows 1-4 are cables
% % Rows 5-10 are bars
% % Columns 1-4 are bottom vertebra nodes
% % Columns 5-8 are top vertebra nodes
% 
% %    1  2  3  4  5  6  7  8  
% C = [0  1  0  0  0 -1  0  0;  %  1, cable 1
%      0  0  1  0  0  0 -1  0;  %  2, ...
%      0  0  0  1  0 -1  0  0;  %  3, ...
%      0  0  0  1  0  0 -1  0;  %  4, cable 4
%      1 -1  0  0  0  0  0  0;  %  5, bar 1
%      1  0 -1  0  0  0  0  0;  %  6, ...
%      1  0  0 -1  0  0  0  0;  %  7, ...
%      0  0  0  0  1 -1  0  0;  %  8, ...
%      0  0  0  0  1  0 -1  0;  %  9, ...
%      0  0  0  0  1  0  0 -1]; % 10, bar 6

% Need to specify number of cables, to split up C.
% Cables per pair of bodies is
sigma = 4;
% and there is one pair of cables per moving body,
s = sigma * b;
% s = 4;

% r follows directly, it's the remainder number of rows.
r = size(C,1) - s;
% ...because C is \in R^{10 x 8}.

% later, we will need number of nodes per body. Could get num col of any of
% the submatrices here
eta = size(Cso, 2);

% number of nodes
n = size(C, 2);

if debugging >= 2
    C
end

% gravitational constant
g = 9.81;

% mass per nodes according to vertebra body
% 4 nodes per body, so
mNode = mBody/4;

% and in vector form,
m = ones(n, 1) * mNode;

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 4, and keep nodes 5-8
% (which is vertebra two.)
% w = [0; 0; 0; 0; 1; 1; 1; 1];

% For a b-body structure that removes the leftmost set of nodes,
% zeros for the first body's nodes, ones for the remainder.

%%%%%%%%%%%%%%%%%%%%%% 2019-07-18: POSSIBLE BUG. Should this be eta?
w = [zeros(sigma, 1); 
     kron(ones(b,1), ones(sigma, 1))];

% Including all nodes:
% w = ones(n,1);

% IMPORTANT! If chosing to remove nodes, must change 'b' also, or else the
% iso optimization problem will FAIL.

%% Trajectory of poses

% Number of poses
% numPts = 10;
% numPts = 30; % was 30 for my dissertation
numPts = 20;

% Initial pose of the spine
% The frames need to be rotated to horizontal from the above
rotation0 = -pi/2;
% Each vertebra is to the right of the previous one
translation0 = [barEndpoint * (3/2);
                0];
            
% A good point to rotate everything from is the origin (center of first
% vertebra)
rotAxisPt = [0; 0];

% let's have the vertebra start horizontal and swing upward
% these are referenced to the MAXIMUM sweep angle: these are what the
% furthest-out vertebra performs.
minSweep = 0;
maxSweep = pi/10;

% get the trajectory:
% NOTE that this is for all bodies, but since we're removing one, have to
% add it back in: b+1.
xi = trajBendMultiBody_2d(translation0, rotation0, minSweep, ...
    maxSweep, rotAxisPt, numPts, b+1);

% % all the positions of each rigid body (expressed as their COM positions
% % and rotation). that's 3 states: [x; z; \gamma] with
% % the angle being an intrinsic rotations.
% 
% % use \xi for the system states.
% xi = zeros(b * 3, 1);
% 
% % for rigid body 1, the fixed one, doesn't move. However, we've defined it
% % in its "vertical" state, so it needs to be rotated by 90 degrees around
% % the y-axis so the saddle cables can align. Let's rotate it clockwise so
% % that node 4 is in +x.
% xi(3) = -pi/2;
% 
% % for rigid body 2, translate out in the +x direction. Translating by one
% % full body length puts the tips exactly in the same plane, so maybe do it
% % to 3/4 of that length.
% % the length of one vert is 2 * bar_endpoint. 
% % x-position is coordinate 1.
% % To make things interesting, let's rotate it a small bit, too.
% %xi(3) = -pi/2 + pi/16;
% 
% xi(4:6) = [     barEndpoint * (3/2);
%                 0;
%                -pi/2 + pi/16];
%             
% if debugging >= 2
%     xi
% end

%% Calculations for the inputs to the core tiso library

% The nodal coordinates (x, z)
% calculate from position trajectory


% coord = getCoord_2d(a, xi, debugging);

if debugging >= 2
    coord
end

% % Split up the coordinates.
% % The tiso core routine expects these to be column vectors, so a
% % transpose is needed also
% x = coord(1, :)';
% y = coord(2, :)';

% We can do these out for each point, along the trajectories.
% initialize the results. There are n nodes (n position vectors.)
x = zeros(n, numPts);
y = zeros(n, numPts);

for i=1:numPts
    % We split up according to rigid body. 
    % Need to adjust again ***TO-DO*** really need to standardize this
    % because there are 3 bodies!
    for k=1:b+1
        % There need to be three states passed in as the second argument to
        % this function. indices are
        xi_ik = xi( (3*(k-1)+1) : 3*k, i);
        coord_ik = getCoord_2d(a, xi_ik, debugging);
        % ...which is of size eta x 1.
        % insert into the x and y vectors,
        % assuming that the bodies are numbered one after the other.
        x( (eta*(k-1)+1) : eta*k, i) = coord_ik(1,:)';
        y( (eta*(k-1)+1) : eta*k, i) = coord_ik(2,:)';
    end
end


% Initialize external forces
px = zeros(n, numPts);
py = zeros(n, numPts);

% Add the gravitational reaction forces for each mass.
% a slight abuse of MATLAB's notation: this is vector addition, no indices
% needed, since py and m are \in R^n.
py = py - m*g;

% The objective function weighting matrix, via the helper(s).
% For 0.5 q_s'*R*q_s
% R = getObj_2norm(s);

%%% Need to calculate this per-step... one R for each timestep.
% For the potential energy objective function,
% need to declare some spring constants
lbin_in_nm = 175.126835;
kappa_i = 4.79 * lbin_in_nm;
% vectorized:
kappa = ones(s,1) * kappa_i;

% then calculate the lengths and get an R for each timestep
lengths = zeros(s, numPts);
% lazy here, don't feel like 3D array indexing
R = {};
for i=1:numPts
    lengths(:,i) = getLen_2d(x(:,i), y(:,i), s, C);
    R{end+1} = getObj_PE(s, lengths(:,i), kappa);
end

% so the objective function weight is
% R = getObj_PE(s, lengths, kappa);

%% Solve the inverse statics problem

% % The struct that contains all mandatory inputs to the ISO problem
% mandInputs.x = x;
% mandInputs.y = y;
% mandInputs.px = px;
% mandInputs.py = py;
% mandInputs.C = C;
% mandInputs.b = b;
% mandInputs.s = s;
% 
% % The optional inputs
% optionalInputs.w = w;
% optionalInputs.R = R;
% optionalInputs.qMin = qMin;
% optionalInputs.debugging = debugging;

% Solve. Rigid body reformulation, two dimensions.
% [fOpt, qOpt, Ab, pb] = rbISO_2d(mandInputs, optionalInputs);

fOpt = zeros(s, numPts);
qOpt = zeros(s, numPts);
% we also need the lengths for calculating u from q*.
% TO-DO: THIS SEEMS REDUNDANT FROM THE ABOVE LENGTHS CALCULATION?
len = zeros(s, numPts);

% finally, the big function call:
for i=1:numPts
    if debugging >= 1
        disp('Iteration:');
        disp(num2str(i));
    end
    % combine parameters into the mandatory and optional arguments to core
    % The struct that contains all mandatory inputs to the ISO problem
    mandInputs_i.x = x(:,i);
    mandInputs_i.y = y(:,i);
    mandInputs_i.px = px(:,i);
    mandInputs_i.py = py(:,i);
    mandInputs_i.C = C;
    mandInputs_i.b = b;
    mandInputs_i.s = s;

    % The optional inputs
    optionalInputs_i.w = w;
    optionalInputs_i.R = R{i};
    optionalInputs_i.qMin = qMin;
    optionalInputs_i.debugging = debugging;
    % including the uMin constraint
    optionalInputs_i.uMin = uMin;
    optionalInputs_i.kappa = kappa;
    
    % quadprog is inside this function.
    [fOpt(:,i), qOpt(:,i), len(:,i), ~, ~] = rbISO_2d(mandInputs_i, ...
        optionalInputs_i);
end

%% Plot the structure, for reference.

% This should make it easier to visualize the results.
% Need to specify "how big" we want the bars to be. A good number is
%radius = 0.02; % meters.
% For the hardware test in fall 2018, a better number, based on the
% hardware itself, is 1.5 inches, which in meters is
%radius = 0.038;
% but that's too big given how the plotter works right now. Choose
% something arbitrarily smaller.
% radius = 0.01;

% Plot.
% plotTensegrity2d(C, x, y, s, radius);

%% Plot the structure, for reference.

% This should make it easier to visualize the results.
% Need to specify "how big" we want the bars to be. A good number is
radius = 0.005; % meters.
% 2019-05-13: The actual vertebra is 4 cm in width, so 
% radius = 0.02;
% Need to modify the plotting script to make this more sensible.

% Plotting: do the first position and the last position.
% Basic command:
% plotTensegrity2d(C, x, y, s, radius);
% We want to get the first and last coordinates.
% Reference configuration state:
% plotTensegrity2d(C, calibratedCoordX, calibratedCoordY, s, radius);
% Initial position:
plotTensegrity2d(C, x(:,1), y(:,1), s, radius);

% Make the axis limits nicer.
ylim([-0.2 0.3]);
xlim([-0.11 0.75]);

% Final position:
plotTensegrity2d(C, x(:,end), y(:,end), s, radius);
ylim([-0.2 0.3]);
xlim([-0.11 0.75]);

% plot the cable tensions too
% legendlabels = {'Horiz. Top', 'Horiz. Bot.', 'Saddle Top', 'Saddle Bot.'};
legendlabels = {'HB', 'HT', 'SB', 'ST'};
plotCableTensions(fOpt, sigma, 80, legendlabels);















%% horizontal_spine_invkin_2d_hardwaretesttrajectory.m
% Copyright Andrew P. Sabelhaus 2018

% This script used the tInvKin libraries to calculate the inverse
% kinematics for a tensegrity spine, defined in 2 dimensions. 
% As the term is used here, 'spine'
% refers to a structure with repeated rigid bodies of the same geometry and
% mass.

% The 'trajectory' script here outputs a sequence of inverse kinematics
% solutions, for use with control.

% Modified from original ..._trajectory script to accomodate for the small
% geometric differences with the hardware test setup.

%% set up the workspace
clear all;
close all;
clc;

% add the core libraries, assumed to be in an adjacent folder.
addpath( genpath('../invkin_core') );
% same for the plotting.
addpath( genpath('../invkin_plotting') );
% and for the trajectories, a subfolder.
addpath( genpath('./invkin_trajectories') );

%% Set up the parameters

% Debugging level.
% 0 = no output except for errors
% 1 = Starting message, results from quadprog
% 2 = Verbose output of status.
debugging = 1;

% If appropriate, output a starting message.
if debugging >= 1
    disp('Starting horizontal spine 2d rigid body inverse kinematics example...');
end

% minimum cable force density
%q_min = 0; % units of N/m, depending on m and g
q_min = 0.5;
%q_min = 0.1;

% Local frames. We're going to use the get_node_coordinates script for both
% "vertebrae", but with the "fixed" one having a different frame to account
% for the slightly different exit points of the cables from the eye hooks.

% For the 2D spine control test, fall 2018, the dimensions of the vertebra
% are as follows:
% bar_endpoint = 4 inches
% that's 4 * 2.54 * 0.01 = 0.1016 meters
bar_endpoint = 0.1016;
% BUT, the center of mass is actually offset, due to some construction
% differences with the nonzero-thickness of the bars.
% The center is at -0.81 inches from the geometric center.
% We'll do the origin of this frame as the geometric center,
% so that puts the center at 
center_offset = 0.81 * 2.54 * 0.01;
center_x = -center_offset;

% The endpoints are still at the bounds of an 8x8 box, which puts the frame
% at the following.
% We're also reorienting the Y-shape so it's in the proper rotation in its
% default configuration.
% 2 = bottom left, 3 = top left, 4 = right.
a_free = [   center_x,              0;
        -bar_endpoint,   -bar_endpoint;
        -bar_endpoint,  bar_endpoint;
        bar_endpoint,              0]';
    
% Now, for the fixed vertebra, we're going to also have four nodes,
% but each of them will only be used as anchor points (and will be thrown
% out of the force balance.)
% It's more difficult to specify the center of that frame based on the
% mechanical designs... let's choose the centerline between the two "tip"
% eyebolts at the origin of that frame. It'll be easier to do the
% translation that way. (there is no rotation for this one.)
% the eyebolts are at 0.75 in from their centerpoint
tip_y = 0.75* 2.54 * 0.01;
% The back parts of the "Y" are at two translations: first, to a reference
% point on the motor plate. Choose this to be the top bolt (it's mated
% aligned in Solidworks.)
back_plate_x = -8.5 * 2.54 * 0.01;
back_plate_y = 3.25 * 2.54 * 0.01;
% Then, on the motor plate, the eye bolt location is translated in X from
% this particular location in the bolt pattern by 1.5in. (No y.)
motor_plate_x = 1.5 * 2.54 * 0.01;
% Therefore, the locations are
back_x = back_plate_x + motor_plate_x;
back_y = back_plate_y;

% Insert these into a full frame.
% Nodes will be numbered from bottom to top.
a_fixed = [     back_x,             -back_y;
                0,                  -tip_y;
                0,                  tip_y;
                back_x,             back_y]';

if debugging >= 2
    a_free
    a_fixed
end

% (e.g., vertebra is in an 8x8 inch box.)

% Mass as measured with a scale on 2018-11-18 is about 500g
m_i = 0.495;

% number of rigid bodies
%b = 2;
% When removing the anchor nodes, it's like removing one of the bodies:
b = 1;

% Configuration matrix for WHOLE STRUCTURE. (in future, pattern out.)

% Full connectivity matrix
% Rows 1-4 are cables
% Rows 5-10 are bars
% Columns 1-4 are bottom tetra nodes
% Columns 5-8 are top tetra nodes
% This is different than the original example, since our nodes for the
% fixed vertebra are different. Same ordering of cables though.
%    1  2  3  4  5  6  7  8  
C = [1  0  0  0  0 -1  0  0;  %  1, cable 1, horizontal bottom
     0  0  0  1  0  0 -1  0;  %  2, cable 2, horizontal top
     0  1  0  0  0 -1  0  0;  %  3, cable 3, saddle bottom
     0  0  1  0  0  0 -1  0;  %  4, cable 4, saddle top
     1 -1  0  0  0  0  0  0;  %  5, bar 1, bottom mounting point connection
     0  1 -1  0  0  0  0  0;  %  6, ...
     0  0  1 -1  0  0  0  0;  %  7, ...
     0  0  0  0  1 -1  0  0;  %  8, ..., free vertebra bar 1
     0  0  0  0  1  0 -1  0;  %  9, ...
     0  0  0  0  1  0  0 -1]; % 10, bar 6

% Need to specify number of cables, to split up C.
s = 4;
% r follows directly, it's the remainder number of rows.
r = size(C,1) - s;
% ...because C is \in R^{10 x 8}.

% number of nodes
n = size(C, 2);

if debugging >= 2
    C
end

% gravitational constant
% ON 2018-12-6: THIS MAY BE OPPOSITE? WHY?? OR NOT??
% TO-DO: FIX THIS SIGN ERROR. Gravity went "up" for earlier results??
g = 9.81;
%g = -9.81;

% mass per vertebra
% let's say each vertebra weighs 0.8 kg. Thta's about 1.7 lbs, which seems
% right to Drew if motors are included.
%m_i = 0.8;
% 2018-11-28: updated above alongside the dimensions for the hardware test
% of fall 2018.

% Note here that I've used m_i as per-body not per-node.
% Probably accidentally changed notation w.r.t. T-CST 2018 paper.
% 4 nodes per body, so
m_node = m_i/4;

% and in vector form,
%m = ones(n, 1) * m_node;

% 2018-12-6: let's see if we just assign all the mass to the center of mass
% node. For the moving vertebra, that's node 1, e.g. 5 overall.
m = zeros(n,1);
m(5) = m_i;
% and we can distribute the masses over the first vertebra - not that it
% matters, since we're throwing out those nodes.
m(1:4) = ones(4,1) * m_node;

% Spring constant for the cables. This isn't used in the optimization for
% force (densities) itself, but for conversion into inputs when saving the
% data.
% Here's how to specify 'the same spring constant for all cables'. In N/m.
% For the hardware test, the Jones Spring Co. #241 has a constant of 1.54
% lb/in, which is 270 N/m.
% One set of new McMaster springs is 4.79 lb/in,
% conversion is 
lbin_in_nm = 175.126835;
kappa_i = 4.79 * lbin_in_nm;
% ...which happens to be
%kappa_i = 270;
kappa = ones(s,1) * kappa_i;

% Example of how to do the 'anchored' analysis.
% Declare a vector w \in R^n, 
% where w(i) == 1 if the node should be 'kept'.
% For this example, want to treat body 1 as the anchored nodes.
% So, we zero-out anchored nodes 1 through 4, and keep nodes 5-8
% (which is vertebra two.)
w = [0; 0; 0; 0; 1; 1; 1; 1];
% Including all nodes:
%w = ones(n,1);

% IMPORTANT! If chosing to remove nodes, must change 'b' also, or else inv
% kin will FAIL.

% We also need to declare which nodes are pinned (with external reaction
% forces) and which are not.
% We're choosing not to have this be the same as w, since there are some
% nodes to "ignore" for w which are not necessarily built in to the ground
% for "pinned". Example, nodes inside the leftmost vertebra, where we're
% deciding to assume that only the tips of its "Y" are supported.

% So, only nodes 1 and 4 are pinned.
pinned = zeros(n,1);
pinned(1) = 1;
pinned(4) = 1;

%% Trajectory of positions

% all the positions of each rigid body (expressed as their COM positions
% and rotation). that's 3 states: [x; z; \gamma] with
% the angle being an intrinsic rotation.

% The trajectory generation function gives back a sequence of states for
% which we require cable inputs.
% Let's do a trajectory that sweeps from 0 to pi/8.
%max_sweep = pi/8;
max_sweep = pi/12;
% with the the vertebra horizontal-sideways,
% The local frame needs to be rotated by
%rotation_0 = -pi/2;
% The frames are at an initial configuration of zero rotation.
rotation_0 = 0;
% and translated outward. The test setup's default position has the free
% vertebra directly aligned with the saddle cables' eyebolt exit points,
% which, since the moving vertebra has its frame origin at the geometric
% center,
% translation_0 = [bar_endpoint; 0];
% For testing, here's a smaller fraction.
translation_0 = [bar_endpoint * (3/4); 0];
% with a large number of points.
num_points = 400;

% For the frames we've chosen, it doesn't make sense to rotate the vertebra
% around the origin: we don't want it to sweep out from the tip of the
% eyebolts, but instead from whatever the "center" of the fixed vertebra
% eyebolt pattern would be.
% That's maybe halfway between the tip and the back, or
rot_axis_pt = [back_x*(2/3); 0];

% get the trajectory:
[xi_all] = trajectory_XZT_bend_2d(translation_0, rotation_0, max_sweep, rot_axis_pt, num_points);

% use \xi for the system states.
%xi = zeros(b * 3, 1);

% for rigid body 1, the fixed one, doesn't move. However, we've defined it
% in its "vertical" state, so it needs to be rotated by 90 degrees around
% the y-axis so the saddle cables can align. Let's rotate it clockwise so
% that node 4 is in +x.
% xi(3) = -pi/2;

% for rigid body 2, translate out in the +x direction. Translating by one
% full body length puts the tips exactly in the same plane, so maybe do it
% to 3/4 of that length.
% the length of one vert is 2 * bar_endpoint. 
% x-position is coordinate 1.
% To make things interesting, let's rotate it a small bit, too.
%xi(3) = -pi/2 + pi/16;

% xi(4:6) = [     bar_endpoint * (3/2);
%                 0;
%                -pi/2 + pi/16];
%             
% if debugging >= 2
%     xi
% end

%% Calculations for the inputs to the core invkin library

% The nodal coordinates (x, z)
% calculate from position trajectory. 
% We can do these out for each point, along the trajectories.
% initialize the results. There are n nodes (n position vectors.)
x = zeros(n, num_points);
y = zeros(n, num_points);

for i=1:num_points
    % We split up according to rigid body. The fixed body is xi_all(1:3,i)
    coordinates_i_fixed = get_node_coordinates_2d(a_fixed, ...
        xi_all(1:3, i), debugging);
    % and for the free vertebra
    coordinates_i_free = get_node_coordinates_2d(a_free, ...
        xi_all(4:6, i), debugging);
    % At this point along the trajectory, get the coordinates
    %coordinates_i = get_node_coordinates_2d(a_free, xi_all(:,i), debugging);
    % ...and split them into coordinate-wise vectors per node, per
    % vertebra.
    x(1:4, i) = coordinates_i_fixed(1,:)';
    x(5:8, i) = coordinates_i_free(1,:)';
    y(1:4, i) = coordinates_i_fixed(2,:)';
    y(5:8, i) = coordinates_i_free(2,:)';
%     x(:,i) = coordinates_i(1,:)';
%     y(:,i) = coordinates_i(2,:)';
end
    
% coordinates = get_node_coordinates_2d(a, xi, debugging);

if debugging >= 2
    %coordinates
    x
    y
end

% Reaction forces can be calculated by this helper, which assumes that only
% gravity is present as an external force.
% Initialize results
px = zeros(n, num_points);
py = zeros(n, num_points);

% Iterate over all points:
for i=1:num_points
    % The i-th index will be columns for all these.
    [px(:,i), py(:,i)] = get_reaction_forces_2d(x(:,i), y(:,i), pinned, ...
        m, g, debugging);
end

% for more details, you can look at commits to the library before Nov. 2018
% where this reaction force/moment balance was written out by hand.

% Since this was just a calculation of the reaction forces, we ALSO need to
% add in the external forces themselves (grav forces) for use as the whole
% inverse kinematics problem.

% Add the gravitational reaction forces for each mass.
% a slight abuse of MATLAB's notation: this is vector addition, no indices
% needed, since py and m are \in R^n.
py = py + -m*g;


%% Solve the inverse kinematics problem

% Solve, over each point.
% Let's use a cell array for Ab and pb, since I don't feel like thinking
% over their sizes right now.
f_opt = zeros(s, num_points);
q_opt = zeros(s, num_points);
% we also need the lengths for calculating u from q*.
lengths = zeros(s, num_points);
Ab = {};
pb = {};

% finally, the big function call:
for i=1:num_points
    if debugging >= 1
        disp('Iteration:');
        disp(num2str(i));
    end
    % quadprog is inside this function.
    [f_opt(:,i), q_opt(:,i), lengths(:,i), Ab_i, pb_i] = invkin_core_2d_rb(x(:,i), ...
        y(:,i), px(:,i), py(:,i), w, C, s, b, q_min, debugging);
    % and insert this Ab and pb.
    Ab{end+1} = Ab_i;
    pb{end+1} = pb_i;
end

% Seems correct, intuitively!
% Cable 1 is horizontal, below.
% Cable 2 is horizontal, above.
% Cable 3 is saddle, below.
% Cable 4 is saddle, above.

% It makes sense that cable 2 force > cable 1 force, for an "upward" force
% in the saggital plane, counteracting gravity on the suspended vertebra.
% it then also could make sense that cable 3 has no force: the
% gravitational force on the suspended vertebra gives a clockwise moment,
% which is also the direction of cable 3.
% Cable 4's force is likely counteracting the gravitational moment of the
% suspending vertebra.

% TO-DO: Confirm, by hand/in simulation, that these forces actually keep
% the vertebra in static equilibrium.

% TO-DO: sign check. Are we applying positive tension forces or negative
% tension forces? See if/what solution pops out with the constraint in the
% opposite direction (< c, not > c.)

% A quick plot of the cable tensions.
figure; 
hold on;
subplot(4,1,1)
hold on;
title('Cable tensions');
plot(f_opt(1,:))
ylabel('1 (N)');
subplot(4,1,2)
plot(f_opt(2,:))
ylabel('2 (N)');
subplot(4,1,3)
plot(f_opt(3,:))
ylabel('3 (N)');
subplot(4,1,4);
plot(f_opt(4,:));
ylabel('4 (N)');

%% Convert the optimal forces into optimal rest lengths.
% u_i = l_i - (F_i/kappa_i)

% save in a vector
u_opt = zeros(s, num_points);
% it's more intuitive to iterate for now. At least, we can iterate over
% cables and not over timesteps.
for k=1:s
    % For cable k, divide the row in f_opt by kappa(k)
    u_opt(k, :) = lengths(k,:) - (f_opt(k,:) ./ kappa(k));
end

% For use with the hardware example, it's easier to instead define a
% control input that's the amount of "stretch" a cable experiences.
% Since this is the definition of F = \kappa * stretch, just obviously the
% amount of spring extension, we can calculate this quantity without even
% using length.
stretch_opt = zeros(s, num_points);
for k=1:s
    % For cable k, divide the row in f_opt by kappa(k)
    % Also, for the hardware test, convert to cm because we're using
    % single-precision floating point numbers.
    stretch_opt(k, :) = (f_opt(k,:) ./ kappa(k)) * 100;
end

%% Plot the structure, for reference.

% This should make it easier to visualize the results.
% Need to specify "how big" we want the bars to be. A good number is
radius = 0.005; % meters.

% Plotting: do the first position and the last position.
% Basic command:
% plot_2d_tensegrity_invkin(C, x, y, s, radius);
% We want to get the first and last coordinates.
% Initial position:
plot_2d_tensegrity_invkin(C, x(:,1), y(:,1), s, radius);
% Final position:
plot_2d_tensegrity_invkin(C, x(:,end), y(:,end), s, radius);

%% Save the data.

% path to store: ***CHANGE THIS PER-USER***
% for now, use the user's home directory.
savefile_path = '~/';

% write the actual data
% we used the rigid body reformulation method here, 
n_or_b = 1;
%save_invkin_results_2d(u_opt, n, r, n_or_b, savefile_path);
% For the hardware test, we want to use "stretch" not rest length.
%save_invkin_results_2d(stretch_opt, n, r, n_or_b, savefile_path);















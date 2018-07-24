% get_node_coordinates_3d.m
% Copyright Andrew P. Sabelhaus

function [coordinates] = get_node_coordinates_3d(local_frame, state, debugging)
% get_node_coordinates calculates the coordinates (positions) of each node
%   for the whole tensegrity structure. For more than 1 rigid body, the
%   'state' vector implies how many times to copy/translate/rotate the
%   local frame.
%   
%   constants taken from inputs:
%   b = number of rigid bodies, from size of state
%   n_i = number of nodes from rigid body, from size of local_frame
%
% Inputs:
%   local_frame = matrix with all points in a local frame for one body.
%       Size 3 x n_i (node coordinates are column vectors.)
%   state = vector, size 6*b where b is the number of rigid bodies.
%       Specifies the translations and rotations for each body.
%
% Outputs:
%   coordinates = matrix, 3 x n_i*b, for all 3 position coordinates for all
%   n_i nodes of b rigid bodies.

% The parameters are
b = size(state,1) / 6;
n_i = size(local_frame, 2);

if debugging
    b
    n_i
end

% Create a matrix in which to place all the coordinates
coordinates = zeros( 3, n_i * b);

% For each rigid body,
for i = 1:b
    % for ease, split off the state for this rigid body
    state_i = state(6*(i-1) + 1 : 6*i);
    if debugging
        state
    end
    % The nodes are at the new xyz coordinates + rotation on the original
    % frame. Rotation matrices are ordered z, y, x
    % debugging
    % A shortcut to the rotation matrices is MATLAB's nice makehgtform
    % command. However, need to truncate off the 4th row/column since the
    % command is used for computer graphics and we don't care about
    % homongeneity.
    rot_z = makehgtform('zrotate', state_i(6));
    rot_z = rot_z(1:3, 1:3);
    rot_y = makehgtform('yrotate', state_i(5));
    rot_y = rot_y(1:3, 1:3);
    rot_x = makehgtform('xrotate', state_i(4));
    rot_x = rot_x(1:3, 1:3);
    
    if debugging
        rot_x
        rot_y
        rot_z
    end

    % apply all the rotations
    rotated_frame = rot_z * rot_y * rot_x * local_frame;
    
    if debugging
        rotated_frame
    end
    
    % Put the translations in matrix form for ease. the frame has each
    % coordinate as a row, so just pattern the translations out as a row
    % vector.
    translations = [ state_i(1) * ones(1, n_i);
                     state_i(2) * ones(1, n_i);
                     state_i(3) * ones(1, n_i)]
        
    % The new frame is then translated also.
    rot_trans_frame = rotated_frame + translations;
    % Finally, insert this into the coordinates matrix.
    % columns indexing is by rigid body, in order. E.g., with 5 points in a
    % local frame and 2 rigid bodies, these should be blocks 1-5 and 6-10.
    coordinates(:, n_i*(i-1) + 1 : n_i*i) = rot_trans_frame;
end

end


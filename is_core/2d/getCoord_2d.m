% getCoord2d.m
% Copyright Andrew P. Sabelhaus

function [coordinates] = getCoord2d(localFrame, state, debugging)
% getCoord2d calculates the coordinates (positions) of each node
%   for the whole tensegrity structure, in two dimensions, given the 
%   states of each body. For more than 1 rigid body, the
%   'state' vector implies how many times to copy/translate/rotate the
%   local frame.
%   
%   constants taken from inputs:
%   b = number of rigid bodies, from size of state
%   n_i = number of nodes from rigid body, from size of local_frame
%
% Inputs:
%   localFrame = matrix with all points in a local frame for one body.
%       Size 2 x n_i (node coordinates are column vectors.)
%
%   state = vector, size 3*b where b is the number of rigid bodies.
%       Specifies the translations and rotations for each body.
%       (3 states are x, z, \gamma.)
%
%   debugging = level of debugging/verbosity from this script. Options:
%       0 or 1 = no output except for errors
%       2 = Verbose output, lots of dimensions and matrices and whatnot.
%
% Outputs:
%   coordinates = matrix, 2 x n_i*b, for all 2 position coordinates for all
%   n_i nodes of b rigid bodies.

% The parameters are
b = size(state,1) / 3;
n_i = size(localFrame, 2);

if debugging >= 2
    b
    n_i
end

% Create a matrix in which to place all the coordinates
coordinates = zeros( 2, n_i * b);

% For each rigid body,
for i = 1:b
    % for ease, split off the state for this rigid body
    state_i = state(3*(i-1) + 1 : 3*i);
    if debugging >= 2
        state_i
    end
    % The nodes are at the new x z coordinates + rotation on the original
    % frame.
    % Rotation matrix for a counterclockwise rotation around the origin is
    gamma = state_i(3);
    rot = [cos(gamma),  -sin(gamma);
           sin(gamma),   cos(gamma)]; 
       
    % apply the rotation
    rotatedFrame = rot * localFrame;
    
    if debugging >= 2
        rotatedFrame
    end
    
    % Put the translations in matrix form for ease. the frame has each
    % coordinate as a row, so just pattern the translations out as a row
    % vector.
    translations = [ state_i(1) * ones(1, n_i);
                     state_i(2) * ones(1, n_i)];
        
    % The new frame is then translated also.
    rotTransFrame = rotatedFrame + translations;
    % Finally, insert this into the coordinates matrix.
    % columns indexing is by rigid body, in order. E.g., with 4 points in a
    % local frame and 2 rigid bodies, these should be blocks 1-4 and 5-8.
    coordinates(:, n_i*(i-1) + 1 : n_i*i) = rotTransFrame;
end

end


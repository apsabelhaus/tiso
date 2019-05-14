% save_invkin_results_2d.m
% Copyright Andrew P. Sabelhaus and BEST Lab, 2018

function save_invkin_results_2d(u_opt, xi, n, r, n_or_b, path)
% save_invkin_results_2d saves the optimal rest lengths for a trajectory of
% inverse kinematics calculations in 2d. 
%
% Inputs:
%   u_opt = an s x num_points matrix of optimal rest lengths, calculated
%       f_opt via the spring constant k. 
%       Or, whatever input would like to be saved (could do stretch).
%   xi = the trajectory of states that were used to calculate u_opt.
%   n = scalar, number of nodes, for a debugging message printed to file
%   s = scalar, number of cables, for a debugging message printed to file
%   r = scalar, number of bars, for a debugging message printed to file
%   n_or_b = {0,1}, where 0 is if the results are calculated via the nodal
%       method and 1 if the results are via the rigid body reformulated
%       method. Again, just for debugging.
%   path = directory into which to save the CSV file.
%
% Outputs:
%   (to MATLAB): no returned data.
%   (to filesystem): a CSV file with the invkin results and some info about
%       the process. Note that the output is transposed in the csv file for
%       ease of use in other programs.

%% some setup

% the number of cables is the number of rows of u_opt.
s = size(u_opt,1);

%% Make the string header for debugging.
s1 = 'Tensegrity Inverse Kinematics Results';
s2 = 'Parameters:';
s3 = 'Dimensionality, Timestamp, Nodes, Cables, Bars, Computation Method';

% Record the current time in a string
start_time_string = datestr(datetime('now'));
% Remove the colons and spaces from this string so that Windows doesn't complain when this repository is cloned
% colons become dashes, spaces become underscores. Regular expressions to the rescue!
start_time_string = regexprep(start_time_string, ':', '-');
start_time_string = regexprep(start_time_string, ' ', '_');

% Need to output rigid body OR nodal method
method = '';
if n_or_b == 0
    method = 'nodal';
elseif n_or_b == 1
    method = 'rigid body';
end

% this is a two-dimensional problem
d = 2;

% second line will be a big list of all the parameters.
s4 = sprintf('%i, %s, %i, %i, %i, %s', d, start_time_string, n, s, r, method);

% describe what the outputs are.
%s5 = 'Optimal rest lengths for each cable starting from cable 1 up to cable s (rows are timestep and columns are cable):';
s5 = 'Inputs: row = timestep and col = cable no.,Col > s = states,States are x y rot (per body)';

%% Write the header

% the filename itself is
filename = strcat('invkin_results_2d_', start_time_string, '.csv');
full_filepath = strcat(path, filename);

% A bit of info for the user
disp('Saving inverse kinematics results to file:');
disp(full_filepath);

% open for writing
fid = fopen(full_filepath, 'w');
% output the headers
fprintf(fid, '%s\n', s1);
fprintf(fid, '%s\n', s2);
fprintf(fid, '%s\n', s3);
fprintf(fid, '%s\n', s4);
fprintf(fid, '%s\n', s5);
% close so we can use csvwrite
fclose(fid);

%% Write the inv kin results themselves

% We need to concatenate the states to the inverse kinematics inputs.
output = [u_opt', xi'];

% csvwrite does this nicely for us.
% actually there's no append option there, so use dlmwrite just with a
% comma as the delimiter.
dlmwrite(full_filepath, output, '-append');

end









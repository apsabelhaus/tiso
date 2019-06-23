% saveISOresult2d.m
% Copyright Andrew P. Sabelhaus and BEST Lab, 2019

function saveISOresult2d(uOpt, xi, n, r, nOrB, path)
% saveISOresult2d saves the output of the inverse statics optimization 
%   problem for a trajectory of calculations in two dimensions. This
%   function will write whatever `uOpt' represents, and thus can be unclear
%   whether the output is rest length of a cable, or amount of stretch in
%   a cable, both similar forms of output.
%   BE CAREFUL when using on hardware, make sure your assumptions match!
%
% Inputs:
%   uOpt = an s x num_points matrix of optimal rest lengths, calculated
%       fOpt via the spring constant k. 
%       Or, whatever input would like to be saved (could do stretch).
%   xi = the trajectory of states that were used to calculate uOpt.
%   n = scalar, number of nodes, for a debugging message printed to file
%   s = scalar, number of cables, for a debugging message printed to file
%   r = scalar, number of bars, for a debugging message printed to file
%   nOrB = {0,1}, where 0 is if the results are calculated via the nodal
%       method and 1 if the results are via the rigid body reformulated
%       method. Again, just for debugging.
%   path = directory into which to save the CSV file.
%
% Outputs:
%   (to MATLAB): no returned data.
%   (to filesystem): a CSV file with the inverse statics otimization results and some info about
%       the process. Note that the output is transposed in the csv file for
%       ease of use in other programs.

%% some setup

% the number of cables is the number of rows of uOpt.
s = size(uOpt,1);

%% Make the string header for debugging.
s1 = 'Tensegrity Inverse Statics Optimization Results';
s2 = 'Parameters:';
s3 = 'Dimensionality, Timestamp, Nodes, Cables, Bars, Computation Method';

% Record the current time in a string
startTimeString = datestr(datetime('now'));
% Remove the colons and spaces from this string so that Windows doesn't complain when this repository is cloned
% colons become dashes, spaces become underscores. Regular expressions to the rescue!
startTimeString = regexprep(startTimeString, ':', '-');
startTimeString = regexprep(startTimeString, ' ', '_');

% Need to output rigid body OR nodal method
method = '';
if nOrB == 0
    method = 'nodal';
elseif nOrB == 1
    method = 'rigid body';
end

% this is a two-dimensional problem
d = 2;

% second line will be a big list of all the parameters.
s4 = sprintf('%i, %s, %i, %i, %i, %s', d, startTimeString, n, s, r, method);

% describe what the outputs are.
%s5 = 'Optimal rest lengths for each cable starting from cable 1 up to cable s (rows are timestep and columns are cable):';
s5 = 'Inputs: row = timestep and col = cable no.,Col > s = states,States are x y rot (per body)';

%% Write the header

% the filename itself is
fileName = strcat('tisoResults2d_', startTimeString, '.csv');
fullFilepath = strcat(path, fileName);

% A bit of info for the user
disp('Saving inverse statics optimization results to file:');
disp(fullFilepath);

% open for writing
fid = fopen(fullFilepath, 'w');
% output the headers
fprintf(fid, '%s\n', s1);
fprintf(fid, '%s\n', s2);
fprintf(fid, '%s\n', s3);
fprintf(fid, '%s\n', s4);
fprintf(fid, '%s\n', s5);
% close so we can use csvwrite
fclose(fid);

%% Write the ISO results themselves

% We need to concatenate the states to the inverse statics inputs.
output = [uOpt', xi'];

% csvwrite does this nicely for us.
% actually there's no append option there, so use dlmwrite just with a
% comma as the delimiter.
dlmwrite(fullFilepath, output, '-append');

end









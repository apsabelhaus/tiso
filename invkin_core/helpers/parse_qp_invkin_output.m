% parse_qp_invkin_output.m
% Copyright Andrew Sabelhaus and Berkeley Emergent Space Tensegrities
% Lab, 2018

function parse_qp_invkin_output(exitflag, output_info)
% parse_qp_invkin_output
%   A parser for the output of quadprog, interpreting it in the context of
%   the inverse kinematics problem for tensegrity robots.
%   This makes it easier than trying to interpret qp's output natively.
%
%   Inputs:
%       exitflag = exit flag returned by quadprog
%       output_info = the 'output' struct returned by quadprog
%
%   Outputs:
%       none. The function does not return anything, but instead writes to
%       the command prompt.
%

disp('Inverse Kinematics Result:');

% Switch on the exit flag so we can interpret results
switch exitflag
    
    case 1
        disp('Result found, valid solution returned.');
        
    case 0
        disp('Maximum iterations exceeded. The solution may not be valid, though is likely close.');
        disp('Manually set the options for quadprog to a higher number of iterations if desired.');
        disp('Number of iterations were:');
        disp(output_info.iterations);
        
    case -2
        disp('No solution found. Either there is no solution to one of the constraints, or there is no solutions that satisfes all.');
        disp('This may mean one of a few things. Possible reasons:');
        disp('1) Equilibrium constraint has no solutions');
        disp('2) Only solutions to equilibrium constraint cause cables to be in compression (tension constraint violation)');
        disp('3) If used/implemented, the input saturation constraint may not be possible.');
        disp('Recommendations:');
        disp('Set verbosity level higher for more information.');
        disp('Set a lower minimum tension, c.');
        disp('Examine output for dimensionality and rank of Ab.');
        
    otherwise
        disp('Quadprog exit flag not recognized, set to a higher debugging level for more information');
        
end

end

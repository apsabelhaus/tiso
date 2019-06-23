% parseISOresult.m
% Copyright Andrew Sabelhaus and Berkeley Emergent Space Tensegrities
% Lab, 2019

function parseISOresult(exitFlag, outputInfo)
% parseISOresult
%   A parser for the output of quadprog, interpreting it in the context of
%   the inverse statics optimization program result.
%   This makes it easier than trying to interpret qp's output natively.
%
%   Inputs:
%       exitFlag = exit flag returned by quadprog
%       outputInfo = the 'output' struct returned by quadprog
%
%   Outputs:
%       none. The function does not return anything, but instead writes to
%       the command prompt.
%

disp('Inverse Statics Optimization Result:');

% Switch on the exit flag so we can interpret results
switch exitFlag
    
    case 1
        disp('Result found, valid solution returned.');
        
    case 0
        disp(' ');
        warning('Maximum iterations exceeded. The solution may not be valid, and is likely invalid.');
        disp(' ');
        disp('This may be caused by many conditions, but in particular, may happen if the equilibrium constraint matrix is rank-deficient but the constraint is inconsistent.');
        disp('For example, if A is tall but has rank less than size(q), but rank([A,b]) > rank(A), in which case the static equilibrium eqns are inconsistent, and no solution exists.')
        disp('Consider removing anchors, if possible.');
        disp(' ');
        disp('Alternatively - may be a numerical issue, so you may try manually setting the options for quadprog to a higher number of iterations if desired.');
        disp('Number of iterations were:');
        disp(outputInfo.iterations);
        
    case -2
        disp(' ');
        disp('No solution found. Either there is no solution to one of the constraints, or there is no solutions that satisfes all.');
        disp('This may mean one of a few things. Possible reasons:');
        disp('1) Equilibrium constraint has no solutions');
        disp('2) Only solutions to equilibrium constraint cause cables to be in compression (tension constraint violation)');
        disp('3) If used/implemented, the input saturation constraint may not be possible.');
        disp(' ');
        disp('Recommendations:');
        disp('Set verbosity level higher for more information.');
        disp('Set a lower minimum tension, c.');
        disp('Examine output for dimensionality and rank of equilibrium matrix.');
        disp('Increase spring constants to give a wider range for input saturation constraint.');
        disp(' ');
        error('Inverse Statics Optimization calculation failed, now exiting.');
        
    otherwise
        disp('Quadprog exit flag not recognized, set to a higher debugging level for more information');
        
end

end

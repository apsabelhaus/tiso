% Copyright Albert Li 2018

% Constructs the moment equation matrix from coordinate vectors
function M = get_M(n,x,y,z)
    % NOTE: This isn't ready for dynamics yet, I'm not subtracting off any
    % reference coordinates to get rotations about COM. This should work
    % fine for statics.
    R = diag(x);
    S = diag(y);
    T = diag(z);
    
    M = [zeros(n) -T       S;
         T        zeros(n) -R;
         -S       R        zeros(n)];
end


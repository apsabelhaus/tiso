% Copyright Albert Li and Andrew Sabelhaus 2018-2019

% Constructs the moment equation matrix in 3D from coordinate vectors
function B = getB_3d(n,x,y,z)
    % NOTE: This isn't ready for dynamics yet, I'm not subtracting off any
    % reference coordinates to get rotations about COM. This should work
    % fine for statics.
    X = diag(x);
    Y = diag(y);
    Z = diag(z);
    
    B = [zeros(n) -Z       Y;
         Z        zeros(n) -X;
         -Y       X        zeros(n)];
end


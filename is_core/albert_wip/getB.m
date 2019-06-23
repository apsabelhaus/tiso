function B = getB(n,x,y,z)
    % R, S, T will have subtracted term in dynamic application
    R = diag(x);
    S = diag(y);
    T = diag(z);
    
    B = [zeros(n) -T       S;
         T        zeros(n) -R;
         -S       R        zeros(n)];
end

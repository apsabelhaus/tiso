function C = getC(Cso,Csi,Crv,b)
    k1 = [eye(b-1) zeros(b-1,1)]; % source diagonals
    k2 = [zeros(b-1,1) eye(b-1)]; % sink diagonals
    
    % cable/rod connectivity matrices
    Cs = kron(k1,Cso) + kron(k2,Csi);
    Cr = kron(eye(b),Crv);
    
    C = [Cs;Cr];
end

function Ab = getAb(A,B,K,H)
    Af = K*A*H'; % forces
    Am = K*B*A*H'; % moments
    
    Ab = [Af;Am];
end

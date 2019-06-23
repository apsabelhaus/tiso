% Copyright Albert Li 2018

% Gets the Ab matrix from constituent intermediates
function Ab = get_Ab(d,b,n,A,M,Hs)
    eta = n/b; % nodes/body
    
    K = kron(eye(d*b),ones(1,eta));
    
    % force/moment equation construction
    Af = K*A*Hs';
    Am = K*M*A*Hs';
    
    Ab = [Af;Am];
end

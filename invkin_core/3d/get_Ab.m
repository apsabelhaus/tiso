% Copyright Albert Li 2018

% Gets the Ab matrix from constituent intermediates
function Ab = get_Ab(d,b,n,A,M,Hs)
    eta = n/b; % nodes/body
    
    bodyCollapseMat = kron(eye(d*b),ones(1,eta));
    
    % force/moment equation construction
    Af = bodyCollapseMat*A*Hs';
    Am = bodyCollapseMat*M*A*Hs';
    
    Ab = [Af;Am];
end

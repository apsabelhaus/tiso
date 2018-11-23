% Copyright Albert Li 2018

% Calculates the node-branch incidence matrix C for the class of spinal
% tensegrity comprising Belka
function C = get_tetrahedral_spine_C_3d(b)
    % Cable Source Matrix
    Cso = [1 0 0 0 0; % TEST CABLES
           1 0 0 0 0; % CENTRAL VERT NODE CONNECTED TO NODES 2/3
        
           0 1 0 0 0;
           0 0 1 0 0;
           0 0 0 1 0;
           0 0 0 0 1;
           
           0 0 0 1 0;
           0 0 0 1 0;
           0 0 0 0 1;
           0 0 0 0 1];
    
    % Cable Sink Matrix
    Csi = [0 -1 0 0 0; % TEST CABLES
           0 0 -1 0 0;
        
           0 -1 0 0 0;
           0 0 -1 0 0;
           0 0 0 -1 0;
           0 0 0 0 -1;
           
           0 -1 0 0 0;
           0 0 -1 0 0;
           0 -1 0 0 0;
           0 0 -1 0 0];
    
    kronmat1 = [eye(b-1) zeros(b-1,1)];
    kronmat2 = [zeros(b-1,1) eye(b-1)];
    Cs = kron(kronmat1,Cso) + kron(kronmat2,Csi);
    
    % Bar Matrix
    Crv = [1 -1 0 0 0;
           1 0 -1 0 0;
           1 0 0 -1 0;
           1 0 0 0 -1];
    Cr = kron(eye(b),Crv);
    
    C = [Cs;Cr];
end


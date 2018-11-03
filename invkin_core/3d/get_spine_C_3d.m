% Copyright Albert Li 2018

% Calculates the node-branch incidence matrix C for the class of spinal
% tensegrity comprising Belka
function C = get_spine_C_3d(b)
    % Vertical Cables. Constructs repeating blocks for this spine
    cvblk1 = [0 1 0 0 0;
               0 0 1 0 0;
               0 0 0 1 0;
               0 0 0 0 1];
    cvblk2 = [0 -1 0 0 0;
               0 0 -1 0 0;
               0 0 0 -1 0;
               0 0 0 0 -1];
    kronmat1 = [eye(b-1) zeros(b-1,1)];
    kronmat2 = [zeros(b-1,1) eye(b-1)];
    Cv = kron(kronmat1,cvblk1) + kron(kronmat2,cvblk2);
    
    % Saddle Cables
    csblk1 = [0 0 0 1 0;
               0 0 0 1 0;
               0 0 0 0 1;
               0 0 0 0 1];
    csblk2 = [0 -1 0 0 0;
               0 0 -1 0 0;
               0 -1 0 0 0;
               0 0 -1 0 0];
    Csad = kron(kronmat1,csblk1) + kron(kronmat2,csblk2);
    
    Cs = [Cv;Csad]; % Constructed Cs
    
    % Bars
    crblk = [1 -1 0 0 0;
             1 0 -1 0 0;
             1 0 0 -1 0;
             1 0 0 0 -1];
    Cr = kron(eye(b),crblk); % Constructed Cr
    
    C = [Cs;Cr];
end


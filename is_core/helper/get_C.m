% Albert H. Li and Andrew P. Sabelhaus 2019

function C = get_C(Cso,Csi,Crb,b)
%% get_C
%
%   getC returns the connectivity matrix for a multi-body tensegrity, by
%   patterning out the source, sink, and body matrices.
%
%   here, \eta = nodes per body,
%         \sigma = cables per pair of bodies
%         \rho = bars per body
%
%   Inputs:
%       Cso = 'source' matrix for cable connections. This should be size
%           sigma x eta.
%
%       Csi = 'sink' matrix for cable connections. Same size as Cso.
%
%       Crb = rod matrix per body. Size rho x eta.
%
%       b = number of bodies to pattern out. Must be > 1.
%
%   Outputs:
%       C = full connectivity matrix for the entire structure
%


k1 = [eye(b-1) zeros(b-1,1)]; % source diagonals
    k2 = [zeros(b-1,1) eye(b-1)]; % sink diagonals
    
    % cable/rod connectivity matrices
    Cs = kron(k1,Cso) + kron(k2,Csi);
    Cr = kron(eye(b),Crb);
    
    C = [Cs;Cr];
end

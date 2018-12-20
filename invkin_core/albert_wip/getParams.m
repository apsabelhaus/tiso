function [Cso,Csi,Crv,eta,sigma,rho,n,s,r,H] = getParams(spine,b)
    Cso = spine.Cso;
    Csi = spine.Csi;
    Crv = spine.Crv;
    eta = spine.eta;
    sigma = spine.sigma;
    rho = spine.rho;
    n = eta*b;
    s = sigma*(b-1);
    r = rho*b;
    H = [eye(s) zeros(s,r)];
end

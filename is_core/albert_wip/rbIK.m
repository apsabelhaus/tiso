function qs = rbIK(mi,oi)
    % extracting necessary mandatory inputs
    Cso = mi.spine.Cso;
    Csi = mi.spine.Csi;
    Crv = mi.spine.Crv;
    eta = mi.eta;
    b = mi.b;
    x = mi.x;
    y = mi.y;
    z = mi.z;
    n = mi.n;
    p = mi.p;
    s = mi.s;
    W = mi.W;
    H = mi.H;
    rbcostfunc = mi.rbcostfunc;

    % extracting optional inputs

    % validating inputs
    validate(mi,oi);
    
    % solving optimal member forces
    C = getC(Cso,Csi,Crv,b);
    An = getAn(C,W,x,y,z);
    K = getK(b,eta);
    B = getB(n,x,y,z);
    Ab = [K*An*H';K*B*An*H'];
    beq = [K*p;K*B*p];
    qs = rbcostfunc(Ab,beq);
    
    % plot
    plot_3d_tensegrity_invkin(C,s,W,x,y,z)
end

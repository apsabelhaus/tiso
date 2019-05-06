function q = nodalIK(mi,oi)
    % extracting necessary mandatory inputs
    Cso = mi.spine.Cso;
    Csi = mi.spine.Csi;
    Crv = mi.spine.Crv;
    b = mi.b;
    x = mi.x;
    y = mi.y;
    z = mi.z;
    s = mi.s;
    W = mi.W;
    ncostfunc = mi.ncostfunc;

    % extracting optional inputs

    % validating inputs
    validate(mi,oi);
    
    % solving optimal member forces
    C = getC(Cso,Csi,Crv,b);
    An = getAn(C,W,x,y,z);
    q = ncostfunc(An);
    
    % plot
    plot_3d_tensegrity_invkin(C,s,W,x,y,z)
end

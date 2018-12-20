function validate(mi,oi)

   if(any(structfun(@isempty,mi)) == 1)
       error('You are missing mandatory inputs!')
   end
   
   % extracting mandatory inputs
   Cso = mi.spine.Cso;
   Csi = mi.spine.Csi;
   Crv = mi.spine.Crv;
   eta = mi.spine.eta;
   n = mi.n;
   d = mi.d;
   b = mi.b;
   p = mi.p;
   U = mi.U;
   V = mi.V;
   W = mi.W;
   
   % intermediates
   ccomp = [size(Cso,2) size(Csi,2) size(Crv,2)];
   [csor,~,csov] = find(Cso);
   [csir,~,csiv] = find(Csi);
   [crvpr,~,crvpv] = find(Crv==1);
   [crvnr,~,crvnv] = find(Crv==-1);
   
   if(~isa(d,'double') || (d~=2 && d~=3))
       error('d must be 2 or 3!')
       
   elseif(~isa(b,'double') || b<1)
       error('b must be a positive number!')
       
   elseif(~all(ccomp(:)==ccomp(1)))
       error('Cso, Csi, and Crv must have the same number of columns!')
       
   elseif(length(csor)~=length(unique(csor)) || ~all(csov(:)==csov(1)))
       error('Each row of Cso must have exactly one "1"!')
       
   elseif(length(csir)~=length(unique(csir)) || ~all(csiv(:)==csiv(1)))
       error('Each row of Csi must have exactly one "-1"!')
       
   elseif(length(crvpr)~=length(unique(crvpr)) || ~all(crvpv(:)==crvpv(1)))
       error('Each row of Crv must have exactly one "1" and one "-1"!')
       
   elseif(length(crvnr)~=length(unique(crvnr)) || ~all(crvnv(:)==crvnv(1)))
       error('Each row of Crv must have exactly one "1" and one "-1"!')
       
   elseif(~isequal(size(U),[6 b]))
       error('U must be 6 by b!')
       
   elseif(~isequal(size(V),[3 eta]))
       error('V must be 3 by eta!')
       
   elseif(~isequal(size(W),[d*n d*n])|| ~isdiag(W) || ... 
          ~isempty(find(W~=0 & W~=1, 1)))
       error('W must be dn by dn, diagonal, and contain only "0" and "1"!')
      
   elseif(~isequal(size(p,1),d*n))
       error('p must be length dn!')
   end
   
end

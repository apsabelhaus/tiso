function An = getAn(C,W,x,y,z)
    An = W*[C'*diag(C*x);
            C'*diag(C*y);
            C'*diag(C*z)];
end

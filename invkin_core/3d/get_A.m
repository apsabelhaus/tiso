function A = get_A(C,x,y,z)
    A = [C'*diag(C*x);
         C'*diag(C*y);
         C'*diag(C*z)];
end

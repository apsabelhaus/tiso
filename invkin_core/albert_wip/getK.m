function K = getK(b,eta)
    % assume 3D - dimension parameter passes in correct values for z
    K = kron(eye(3*b),ones(1,eta));
end

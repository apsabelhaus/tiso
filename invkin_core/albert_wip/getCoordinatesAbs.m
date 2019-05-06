function [x,y,z] = getCoordinatesAbs(n,U,V)
    % Input U has ABSOLUTE rotations and translations

    b = size(U,2); % number of bodies
    eta = size(V,2); % nodes per body
    
    coordinates = zeros(3,n);
    
    for body = 1:b
        u = U(:,body); % body's inertial rotation and translation vector
        
        % absolute translations and rotations of body coordinate frame
        xt = u(1);
        yt = u(2);
        zt = u(3);
        xr = u(4);
        yr = u(5);
        zr = u(6);
        
        % retrieving rotation and translation matrices
        R = makehgtform('zrotate',zr,'yrotate',yr,'xrotate',xr);
        R = R(1:3,1:3);
        T = [xt;yt;zt];
        
        % transformation of body nodes
        vNode = R*V + kron(ones(1,eta),T);
        coordinates(:,(eta*(body-1)+1):(eta*body)) = vNode;
    end
    
    x = coordinates(1,:)';
    y = coordinates(2,:)';
    z = coordinates(3,:)';
    
end

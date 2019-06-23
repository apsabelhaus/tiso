function [x,y,z] = getCoordinatesRel(n,U,V)
    % Input U has rotations and translations RELATIVE to previous vertebra

    b = size(U,2); % number of bodies
    eta = size(V,2); % nodes per body
    
    coordinates = zeros(3,n);
    
    vPrev = V;
    
    for body = 1:b
        u = U(:,body); % body's relative rotation and translation vector
        
        % relative translations and rotations of body coordinate frame
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
        vNode = R*vPrev + kron(ones(1,eta),T);
        coordinates(:,(eta*(body-1)+1):(eta*body)) = vNode;
        vPrev = vNode;
    end
    
    x = coordinates(1,:)';
    y = coordinates(2,:)';
    z = coordinates(3,:)';
    
end

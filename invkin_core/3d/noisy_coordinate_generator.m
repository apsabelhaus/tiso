function coordinates = noisy_coordinate_generator(a,b,bar_endpoint,distnoise,angnoise,debugging)
    % Initializing orientation vector
    xi = zeros(b*6,1);
    xi(4:6) = [0; pi/2; 0];
    
    for body = 2:b
        % Noise
        xn7p = normrnd(0,distnoise);
        yn7p = normrnd(0,distnoise);
        zn7p = normrnd(0,distnoise);

        xn7a = normrnd(0,angnoise);
        yn7a = normrnd(0,angnoise);
        zn7a = normrnd(0,angnoise);
        
        % Adding noise and generating spine positions
        xi(6*(body-1)+1:6*(body-1)+6) = [2*(body-1)*bar_endpoint + xn7p;
                             0 + yn7p;
                             0 + zn7p;
                             0 + xn7a;
                             pi/2 + yn7a;
                             0 + zn7a];
    end
    
    coordinates = get_node_coordinates_3d(a, xi, debugging);
end


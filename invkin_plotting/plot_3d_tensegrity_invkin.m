function plot = plot_3d_tensegrity_invkin(C,s,W,x,y,z)
    % Retrieving Parameters
    [m,n] = size(C);
    
    % Start of plot
    plot = figure;
    view(3);
    hold on;
    grid on;
    axis equal;
    
    % Plotting Nodes
    for i = 1:n
        % Logic for plotting anchors vs no anchors
        if W(i,i)==0
            plot3(x(i),y(i),z(i),'^','MarkerSize',5,'Color',[0,1,0]);
            plot3(x(i),y(i),z(i),'.','MarkerSize',5,'Color',[0,0,0]);
        else
            plot3(x(i),y(i),z(i),'.','MarkerSize',15,'Color',[0,0,0]);
        end
    end
    
    % Plotting Cables
    for i = 1:s
        Ci = C(i,:);
        indices = find(Ci);
        i1 = indices(1);
        i2 = indices(2);
        xi = [x(i1) x(i2)];
        yi = [y(i1) y(i2)];
        zi = [z(i1) z(i2)];
        plot3(xi,yi,zi,'Color',[1,0,0]);
    end
    
    % Plotting Bars
    for i = s+1:m
        Ci = C(i,:);
        indices = find(Ci);
        i1 = indices(1);
        i2 = indices(2);
        xi = [x(i1) x(i2)];
        yi = [y(i1) y(i2)];
        zi = [z(i1) z(i2)];
        plot3(xi,yi,zi,'Color',[0,0,0],'LineWidth',2);
    end
    
end

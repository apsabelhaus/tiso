%% Demonstration of infeasibility using Nodal Method
% -------------------------------------------

clc
clear

% Underlying graph + physical parameters
C = [0  1  0  0  0 -1  0  0;  %  1, cable 1
     0  0  1  0  0  0 -1  0;  %  2, ...
     0  0  0  1  0 -1  0  0;  %  3, ...
     0  0  0  1  0  0 -1  0;  %  4, cable 4
     1 -1  0  0  0  0  0  0;  %  5, bar 1
     1  0 -1  0  0  0  0  0;  %  6, ...
     1  0  0 -1  0  0  0  0;  %  7, ...
     0  0  0  0  1 -1  0  0;  %  8, ...
     0  0  0  0  1  0 -1  0;  %  9, ...
     0  0  0  0  1  0  0 -1]; % 10, bar 6
g = 9.81;
c = 0.0001; % min tension
Hhat = [eye(4), zeros(4,6)];
options = optimoptions('quadprog','MaxIterations',1000);

%% Vertical Case
x = [0; 1; -1; 0; 0; 1; -1; 0];
y = [0; -1; -1; 1; 1.5; 0.5; 0.5; 2.5];

% Pin at node 2, roller at node 3, Fx2 = 0 by inspection. Use symmetry.
Ry3 = 4*g;
Ry2 = 4*g;

% Set up constraint matrices
An = [C'*diag(C*x); C'*diag(C*y)];
p = [zeros(8,1); -g.*ones(8,1)];

p(10) = Ry2 - g;
p(11) = Ry3 - g;

Hs = eye(4);
fs = zeros(4,1);
Kf = kron(eye(4),ones(1,4)); % 2 force dimensions
Km = kron(eye(2),ones(1,4)); % 1 moment dimension
B = [diag(-y) diag(x)];
Af = Kf*An*Hhat';
Am = Km*B*An*Hhat';
Ab = [Af;Am];
bineq_m = [Kf*p;Km*B*p];

% Ineq constraints. Ai has a negative applied because the form that
% quadprog takes is Ai*q <= bi, and we want each q >= c, so -q <= -c
Ain = -Hhat;
bi = -ones(4,1) * c;

Aib = -eye(4);

% Optimization Problem
H = 2*eye(10); % so we get q'q instead of 0.5*q'q, doesn't really matter
f = zeros(10,1);
q1 = quadprog(H, f, [], [], An, p,[],[],[], options); % no tension constraint - feas
q2 = quadprog(H, f, Ain, bi, An, p,[],[],[], options); % tension constraint - infeas

qs1 = quadprog(Hs, fs, Aib, bi, Ab, bineq_m, [], [], [], options); % cable forces - feas

plot_2d_tensegrity_invkin(C, x, y, 4, .02);

%% Horizontal Case
x = [0;-1;-1;1;1.5;0.5;0.5;2.5];
y = [0;-1;1;0;0;-1;1;0];

% Structure is horizontal, pin at node 2, roller at node 8
% Parameters:
% m = 1 for all nodes
% Solve the statics problem, Rx2 = 0 by inspection
Ry8 = (8.5/3.5) * g;
Ry2 = (8 - (8.5/3.5)) * g;

% Set up constraint matrices
An = [C'*diag(C*x); C'*diag(C*y)];
p = [zeros(8,1); -g.*ones(8,1)];
p(10) = Ry2 - g;
p(16) = Ry8 - g;

bineq_m = [Kf*p;Km*B*p];

Ai = -Hhat;
bi = -ones(4,1) * c;

% Solve the optimization problem
H = 2*eye(10);
f = zeros(10,1);
q3 = quadprog(H, f, [], [], An, p,[],[],[],options); % no tension constraint - infeas
q4 = quadprog(H, f, Ai, bi, An, p,[],[],[],options); % tension constraint - infeas

qs2 = quadprog(Hs, fs, Aib, bi, Ab, bineq_m, [], [], [], options); % cable forces - feas

plot_2d_tensegrity_invkin(C, x, y, 4, .02);

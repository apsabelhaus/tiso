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
c = 1; % min tension
Hhat = [eye(4), zeros(4,6)];
options = optimoptions('quadprog','MaxIterations',100000);

% Removing anchored points
w = [0;0;0;0;1;1;1;1];
W = diag(w);
Wf = kron(eye(2), W);

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
Af = Kf*Wf*An*Hhat';
Am = Km*B*Wf*An*Hhat';
Ab = [Af;Am];
bineq_m = [Kf*Wf*p;Km*B*Wf*p];

% Ineq constraints. Ai has a negative applied because the form that
% quadprog takes is Ai*q <= bi, and we want each q >= c, so -q <= -c
Ain = -Hhat;
bi = -ones(4,1) * c;

Aib = -eye(4);

% Optimization Problem
H = 2*eye(10); % so we get q'q instead of 0.5*q'q, doesn't really matter
f = zeros(10,1);
q1 = quadprog(H, f, [], [], Wf*An, Wf*p,[],[],[], options); % no tension constraint - feas
q2 = quadprog(H, f, Ain, bi, Wf*An, Wf*p,[],[],[], options); % tension constraint - infeas

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

bineq_m = [Kf*Wf*p;Km*B*Wf*p];

Ai = -Hhat;
bi = -ones(4,1) * c;

% Solve the optimization problem
H = 2*eye(10);
f = zeros(10,1);
q3 = quadprog(H, f, [], [], Wf*An, Wf*p,[],[],[],options); % no tension constraint - infeas @ presolve
q4 = quadprog(H, f, Ai, bi, Wf*An, Wf*p,[],[],[],options); % tension constraint - infeas @ presolve

qs2 = quadprog(Hs, fs, Aib, bi, Ab, bineq_m, [], [], [], options); % cable forces - feas

plot_2d_tensegrity_invkin(C, x, y, 4, .02);

%% 3D 5-Vertebra Case with Anchor Removal
% We should confirm the math here all makes sense. This demonstrates what
% we want, but I don't really know how to explain it.
% (1) Given enough iterations, both the nodal methods fail to converge to a
% solution. THIS NEEDS AN EXPLANATION
% (2) The rigid body method almost immediately produces a solution. My cost
% function combines both the infinity norm of the cable tensions and the
% two-norm. When comparing the values of the infeasible nodal solutions and
% the rigid body method, we can see that the rigid body method produces
% more optimal solutions.

b = 5;
s = 8*(b-1);
r = 4*b;
n = 5*b;
eta = n/b;
d = 3;
c = 1; % min tension
bar_endpoint = 1;
g = 9.81;
debugging = 0;
options = optimoptions('quadprog','MaxIterations',100000);

C = get_tetrahedral_spine_C_3d(b);

w = ones(eta*b,1);
w(1:eta) = zeros(eta,1);
w(end-eta+1:end) = zeros(eta,1);
W = diag(w);
Wf = kron(eye(3),W);
p = [zeros(eta*b,1);zeros(eta*b,1);-g*ones(eta*b,1)]; % don't need reactions b/c of anchor removal
[Hs,~] = get_H(s,r);

a = [   0,              0,              0;
        bar_endpoint,     0,              -bar_endpoint;
        -bar_endpoint,    0,              -bar_endpoint;
        0,              bar_endpoint,     bar_endpoint;
        0,              -bar_endpoint,    bar_endpoint]';
    
% Make points noiser than they are in AbConstructionExample
coordinates = noisy_coordinate_generator(a,b,bar_endpoint,bar_endpoint/5,5*pi/180,debugging);

x = coordinates(1,:)';
y = coordinates(2,:)';
z = coordinates(3,:)';

An = get_A(C,x,y,z);
B = get_M(n,x,y,z);
K = kron(eye(d*b),ones(1,eta));

% Nodal Optimization Setup
H = 2*eye(s+r);
f = zeros(s+r,1);
Ai = -Hs;
bi = -ones(s,1) * c;

% Rigid Body Optimization Setup
Aib = [-eye(s);eye(s)];
Af = K*Wf*An*Hs';
Am = K*B*Wf*An*Hs';
Ab = [Af;Am];
beq_m = [K*Wf*p;K*B*Wf*p];
fs = zeros(s,1);
Hrb = eye(s+1);
Hrb(s+1:s+1) = 0;

[q5, ~, exitflag5, ~] = quadprog(H, f, [], [], Wf*An, Wf*p,[],[],[],options); % nodal, no tension constraints
[q6, ~, exitflag6, ~] = quadprog(H, f, Ai, bi, Wf*An, Wf*p,[],[],[],options); % nodal, tension constraints
[qs3, ~, exitflags3, ~] = quadprog(Hrb, [fs;0], [Aib [zeros(s,1);-ones(s,1)]], [bi;zeros(s,1)], [Ab zeros(2*d*b,1)], beq_m); % rigid body

plot_3d_tensegrity_invkin(C,s,w,x,y,z);

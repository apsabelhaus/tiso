%% Workspace Setup
clear all;
close all;
clc;

% Add the core libraries, assumed to be in an adjacent folder.
addpath(genpath('../invkin_core'));

debugging = 0;

%% Problem Setup
% Parameters
s = 8; % Cables
r = 8; % Bars
n = 10; % Nodes
b = 2; % Bodies
eta = n/b; % Nodes/Body
d = 3; % Dimension (Cartesian space)

% Description: https://arxiv.org/pdf/1806.08868.pdf

% Two vertebrae out of that paper are selected. The nodes corresponding to
% the cores of the vertebrae are 1 and 6, 1-5 are the nodes on the bottom
% vertebra that is fixed. 5-10 are the other vertebra on top. Within the
% groups of 5, nodes 2 and 3 represent the two default bottom oriented
% nodes (in practice, these will be indistinguishable due to symmetry).

% The 4 vertical cables take the first four branch indices, branch 1
% connects nodes 2 and 7, branch 2 connects nodes 3 and 8, etc.

% The 4 saddle cables take the next four branch indices, which connect the
% top two nodes on a bottom vertebra with the bottom two nodes on a top
% vertebra. Here, node 4 is connected to nodes 7 and 8. Node 5 is also
% connected to nodes 7 and 8.

% The 4 bars in the bottom vertebra take the next four indices

% We now construct C and the corresponding H matrix.
Cs = [0 1 0 0 0 -1 0 0 0 0
      0 0 1 0 0 0 -1 0 0 0;
      0 0 0 1 0 0 0 -1 0 0;
      0 0 0 0 1 0 0 0 0 -1;
      0 0 0 1 0 0 -1 0 0 0;
      0 0 0 1 0 0 0 -1 0 0;
      0 0 0 0 1 0 -1 0 0 0;
      0 0 0 0 1 0 0 -1 0 0];
Cr = [1 -1 0 0 0 0 0 0 0 0;
      1 0 -1 0 0 0 0 0 0 0;
      1 0 0 -1 0 0 0 0 0 0;
      1 0 0 0 -1 0 0 0 0 0;
      0 0 0 0 0 1 -1 0 0 0;
      0 0 0 0 0 1 0 -1 0 0;
      0 0 0 0 0 1 0 0 -1 0;
      0 0 0 0 0 1 0 0 0 -1];
C = [Cs;Cr];
Hs = [eye(s) zeros(s,r)];

%% Case 1: Drew's Example Reference Setup

% Map to spatial representation using Drew's setup. a is R', the transpose
% of the representation matrix in graph theory literature. Also specify
% pinned nodes and spec gravity.

bar_endpoint = 0.5; % m
a = [   0,              0,              0;
        bar_endpoint,     0,              -bar_endpoint;
        -bar_endpoint,    0,              -bar_endpoint;
        0,              bar_endpoint,     bar_endpoint;
        0,              -bar_endpoint,    bar_endpoint]';
pinned = zeros(n, 1);
pinned(1:3) = 1; % pin bottom 3
g = 9.81; % m/s^2
m_b = 0.8; % mass of one vertebra
m = m_b/eta * ones(n, 1); % uniform mass distribution, kinda funky

% Set up positions using a xi vector
xi = zeros(b * 6, 1);
xi(4:6) = [0; pi/2; 0];
xi(7:12) = [bar_endpoint * (3/4);
            0;
            0;
            0;
            pi/2;
            0];
coordinates = get_node_coordinates_3d(a, xi, debugging);
x = coordinates(1,:)';
y = coordinates(2,:)';
z = coordinates(3,:)';

% Calculating external force vector
[px, py, pz] = get_reaction_forces_3d(coordinates, pinned, m, g, debugging);
px = -px; py = -py; pz = -pz;
for i=1:n
    pz(i) = pz(i) - m(i)*g;
end
p = [px;py;pz];

% ***Force and Moment Expressions***

% Finding Af and pf
Ao = [C'*diag(C*x);
      C'*diag(C*y);
      C'*diag(C*z)]; % intermediate matrix
Af = kron(eye(d*b),ones(1,eta))*Ao*Hs';
pf = kron(eye(d*b),ones(1,eta))*p;

% Finding Am and Pm

% Use these definitions because we are considering the static case. We
% don't subtract a reference coordinate from every value, which would
% probably coincide with the center of mass.
R = diag(x);
S = diag(y);
T = diag(z);

M = [zeros(n) -T         S;
     T          zeros(n) -R;
     -S         R          zeros(n)];
 
Am = kron(eye(d*b),ones(1,eta))*M*Ao*Hs';
pm = kron(eye(d*b),ones(1,eta))*M*p;

% Combine
Ab = [Af;Am];
pb = [pf;pm];

% ***Analysis***
rankAb = rank(Ab);
dimAb = size(Ab);

% Final print
fprintf('Ab is: \n');
disp(Ab);
fprintf('pb is: \n');
disp(pb);
fprintf(['\nAb has dimensions ' num2str(dimAb(1)) 'x' num2str(dimAb(2))...
      ' and rank ' num2str(rankAb) '.\n']);
fprintf('\nRREF of Ab is \n');
disp(rref(Ab));

% We can see that the rows are pairwise linearly dependent, so this
% structure will have rank be halved. We know that column rank = row rank,
% so we only need to pick one perspective. If column rank drops, we know
% that we have cables that are BOTH force and moment redundant. If row
% rank drops, that implies that forces/moments acting on one body can be
% described entirely as combinations of forces/moments acting on another
% body. For a two-body system, this is fairly intuitive. We're taking the
% moments about a universal origin for a static system, and all cable
% forces acting in one direction on body one must act oppositely in body
% two. This explains the pairwise linear dependence of the rows.

%% Case 2: Slightly Perturbed 2-Vertebra Example
% Now let's run the exact same graph but with a slightly perturbed
% orientation and perform the same rank analysis.

% Noise with std = pi/180 = 6 sigma range of 3 degrees
xn = normrnd(0,pi/180);
yn = normrnd(0,pi/180);
zn = normrnd(0,pi/180);

xi2 = zeros(b * 6, 1);
xi2(4:6) = [0; pi/2; 0];
xi2(7:12) = [bar_endpoint * (3/4);
            0;
            0;
            0 + xn;
            pi/2 + yn;
            0 + zn];
coordinates2 = get_node_coordinates_3d(a, xi2, debugging);
x2 = coordinates2(1,:)';
y2 = coordinates2(2,:)';
z2 = coordinates2(3,:)';
[px2, py2, pz2] = get_reaction_forces_3d(coordinates2, pinned, m, g, debugging);
px2 = -px2; py2 = -py2; pz2 = -pz2;
for i=1:n
    pz2(i) = pz2(i) - m(i)*g;
end
p2 = [px2;py2;pz2];

Ao2 = [C'*diag(C*x2);
       C'*diag(C*y2);
       C'*diag(C*z2)]; % intermediate matrix
Af2 = kron(eye(d*b),ones(1,eta))*Ao2*Hs';
pf2 = kron(eye(d*b),ones(1,eta))*p;

R2 = diag(x2);
S2 = diag(y2);
T2 = diag(z2);

M2 = [zeros(n) -T2         S2;
     T2          zeros(n) -R2;
     -S2         R2          zeros(n)];
 
Am2 = kron(eye(d*b),ones(1,eta))*M2*Ao2*Hs';
pm2 = kron(eye(d*b),ones(1,eta))*M2*p2;

Ab2 = [Af2;Am2];
pb2 = [pf2;pm2];

rankAb2 = rank(Ab2);
dimAb2 = size(Ab2);

fprintf('Ab2 is: \n');
disp(Ab2);
fprintf('pb2 is: \n');
disp(pb2);
fprintf(['\nAb2 has dimensions ' num2str(dimAb2(1)) 'x' num2str(dimAb2(2))...
      ' and rank ' num2str(rankAb2) '.\n']);
fprintf('\nRREF of Ab2 is \n');
disp(rref(Ab2));

% This produces the same result, matching our intuition

%% Case 3: 2-Vertebra Example With Node Alignment

% Now we modify the xi vector to align nodes 4+5 with 7+8 (vertebrae tips)

xi3 = zeros(b * 6, 1);
xi3(4:6) = [0; pi/2; 0];
xi3(7:12) = [2*bar_endpoint;
            0;
            0;
            0;
            pi/2;
            0];
coordinates3 = get_node_coordinates_3d(a, xi3, debugging);
x3 = coordinates3(1,:)';
y3 = coordinates3(2,:)';
z3 = coordinates3(3,:)';
[px3, py3, pz3] = get_reaction_forces_3d(coordinates3, pinned, m, g, debugging);
px3 = -px3; py3 = -py3; pz3 = -pz3;
for i=1:n
    pz3(i) = pz3(i) - m(i)*g;
end
p3 = [px3;py3;pz3];

Ao3 = [C'*diag(C*x3);
       C'*diag(C*y3);
       C'*diag(C*z3)]; % intermediate matrix
Af3 = kron(eye(d*b),ones(1,eta))*Ao3*Hs';
pf3 = kron(eye(d*b),ones(1,eta))*p;

R3 = diag(x3);
S3 = diag(y3);
T3 = diag(z3);

M3 = [zeros(n) -T3         S3;
     T3          zeros(n) -R3;
     -S3         R3          zeros(n)];
 
Am3 = kron(eye(d*b),ones(1,eta))*M3*Ao3*Hs';
pm3 = kron(eye(d*b),ones(1,eta))*M3*p3;

Ab3 = [Af3;Am3];
pb3 = [pf3;pm3];

rankAb3 = rank(Ab3);
dimAb3 = size(Ab3);

fprintf('Ab3 is: \n');
disp(Ab3);
fprintf('pb3 is: \n');
disp(pb3);
fprintf(['\nAb3 has dimensions ' num2str(dimAb3(1)) 'x' num2str(dimAb3(2))...
      ' and rank ' num2str(rankAb3) '.\n']);
fprintf('\nRREF of Ab3 is \n');
disp(rref(Ab3));

% Same result. I highly suspect that the rank of Ab is dependent entirely
% on the rigid geometry of the vertebrae and the choice of nodal routing
% points. In other words, the rank condition is entirely geometric. I need
% to prove this.

%% Case 4: Completely Interfering Vertebrae
xi4 = zeros(b * 6, 1);
xi4(4:6) = [0; pi/2; 0];
xi4(7:12) = [0;
             0;
             0;
             0;
             pi/2;
             0];
coordinates4 = get_node_coordinates_3d(a, xi4, debugging);
x4 = coordinates4(1,:)';
y4 = coordinates4(2,:)';
z4 = coordinates4(3,:)';
[px4, py4, pz4] = get_reaction_forces_3d(coordinates4, pinned, m, g, debugging);
px4 = -px4; py4 = -py4; pz4 = -pz4;
for i=1:n
    pz4(i) = pz4(i) - m(i)*g;
end
p4 = [px4;py4;pz4];

Ao4 = [C'*diag(C*x4);
       C'*diag(C*y4);
       C'*diag(C*z4)]; % intermediate matrix
Af4 = kron(eye(d*b),ones(1,eta))*Ao4*Hs';
pf4 = kron(eye(d*b),ones(1,eta))*p;

R4 = diag(x4);
S4 = diag(y4);
T4 = diag(z4);

M4 = [zeros(n) -T4         S4;
     T4          zeros(n) -R4;
     -S4         R4          zeros(n)];
 
Am4 = kron(eye(d*b),ones(1,eta))*M4*Ao4*Hs';
pm4 = kron(eye(d*b),ones(1,eta))*M4*p4;

Ab4 = [Af4;Am4];
pb4 = [pf4;pm4];

rankAb4 = rank(Ab4);
dimAb4 = size(Ab4);

fprintf('Ab4 is: \n');
disp(Ab4);
fprintf('pb4 is: \n');
disp(pb4);
fprintf(['\nAb4 has dimensions ' num2str(dimAb4(1)) 'x' num2str(dimAb4(2))...
      ' and rank ' num2str(rankAb4) '.\n']);
fprintf('\nRREF of Ab4 is \n');
disp(rref(Ab4));

% Unsurprisingly, this matrix drops rank. It's now rank 5. My question is:
% why doesn't it drop more rank? I need to think about the physical
% intuition behind this a little. From looking at the Ab matrix, we've
% seemingly coupled the Z moments and the y forces together, which is why
% mathematically the rank drops by one. I'm trying to further interpret
% this statement.

%% Extra: Connectivity Analysis and Graph Coloring
% Rigid body analysis from the Cr submatrix
Lr = Cr'*Cr;
[Vr, Dr] = eig(Lr);

% You will see that the two lowest eigenvalues of Lr are 0. Their
% corresponding eigenvectors represent rigid bodies in the tensegrity, and
% the nonzero values correspond to node indices belonging to that body. You
% can see this is nodes 1-5 and 6-10 as expected by our manual labeling.

% Same analysis, this time jumbling the node labels by mixing columns
jumbleCr = Cr(:,randperm(size(Cr, 2)));
while isequal(jumbleCr,Cr)
    jumbleCr = Cr(:,randperm(size(Cr, 2)));
end
jLr = jumbleCr'*jumbleCr;
[Vjr, Djr] = eig(jLr);

% You see the eigenvectors corresponding to the 0 eigenspace are different
% now because our node labeling is different. This implies given some
% labeling of which members are cables vs bars, we can always construct the
% form C = [Cs;Cr] from any arbitrarily constructed C_arb.

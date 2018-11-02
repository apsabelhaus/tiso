%% Exposition and Motivation

% The static equilibrium matrix will always admit solutions for our
% optimization problem if it is wider than it is tall and tall matrices
% will still admit solutions if they are rank deficient.

% Certain topologies are thus intractable from standard nodal analysis
% methods. For example, in the case of our spinal structure, each vertebra
% consists of 5 nodes, which adds 15 rows in Cartesian space. Each
% additional vertebra adds 4 bars + 8 cables = 12 members.

% Clearly, the rate at which the dimension of nd space increases with
% bodies is greater than the rate at which member-space increases, so for
% spinal tensegrities with this type of structure, we can never certainly
% admit solutions.

% Furthermore, we can also discount the effect of just adding bars into the
% system. While this is useful for static analysis, this is not useful for
% dynamic robot systems because they do not add degrees of actuation. Thus,
% the problem will extend into dynamic analysis infinitely for any spinal
% robots using this vertebra structure. Otherwise, we could just trivially
% transform the system into a solvable one by arbitrarily picking points in
% the middle of bars and designating them as nodes.

% The goal of this analysis is to determine the conditions in which we can
% design robots with different formulations of equilibrium matrices such
% that we can abuse rank deficiency and thus produce optimal solutions.

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
Cs = [0 1 0 0 0 0 -1 0 0 0; % vertical cables
      0 0 1 0 0 0 0 -1 0 0;
      0 0 0 1 0 0 0 0 -1 0;
      0 0 0 0 1 0 0 0 0 -1;
      
      0 0 0 1 0 0 -1 0 0 0; % saddle cables
      0 0 0 1 0 0 0 -1 0 0;
      0 0 0 0 1 0 -1 0 0 0;
      0 0 0 0 1 0 0 -1 0 0];

Cr = [1 -1 0 0 0 0 0 0 0 0; % body 1
      1 0 -1 0 0 0 0 0 0 0;
      1 0 0 -1 0 0 0 0 0 0;
      1 0 0 0 -1 0 0 0 0 0;
      
      0 0 0 0 0 1 -1 0 0 0; % body 2
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

% Examining the RREF: the 4th vertical cable and 8th saddle cable are
% force/moment redundant.

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

% This produces the same rank result, matching our intuition. However, note
% that the form of the RREF for this case is completely different! With
% some perturbations, all 4 vertical cables are independent, and 2 saddle
% cables are independent.

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

% Same result as for Ab1. This isn't super surprising - there's nothing
% inherently special about lining up the nodes we chose.

%% Case 4: Completely Interfering Vertebrae

% The two vertebrae will be stacked on top of each other.

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

% Unsurprisingly, this matrix drops rank (it is now rank 4). Based on the
% RREF, it looks like we're at a singularity where only the saddle cables
% are contributing to the cable forces required for this configuration. We
% can't rely on rank reduction in these types of configurations, though,
% since this isn't physically possible.

%% Case 5: Positional AND Angular Perturbations

% Same as Ab3 but now we add noise to the position of vertebra 2.

xn5p = normrnd(0,bar_endpoint/100);
yn5p = normrnd(0,bar_endpoint/100);
zn5p = normrnd(0,bar_endpoint/100);

xn5a = normrnd(0,pi/180);
yn5a = normrnd(0,pi/180);
zn5a = normrnd(0,pi/180);

xi5 = zeros(b * 6, 1);
xi5(4:6) = [0; pi/2; 0];
xi5(7:12) = [2*bar_endpoint + xn5p;
            0 + yn5p;
            0 + zn5p;
            0 + xn5a;
            pi/2 + yn5a;
            0 + zn5a];
coordinates5 = get_node_coordinates_3d(a, xi5, debugging);
x5 = coordinates5(1,:)';
y5 = coordinates5(2,:)';
z5 = coordinates5(3,:)';
[px5, py5, pz5] = get_reaction_forces_3d(coordinates5, pinned, m, g, debugging);
px5 = -px5; py5 = -py5; pz5 = -pz5;
for i=1:n
    pz5(i) = pz5(i) - m(i)*g;
end
p5 = [px5;py5;pz5];

Ao5 = [C'*diag(C*x5);
       C'*diag(C*y5);
       C'*diag(C*z5)]; % intermediate matrix
Af5 = kron(eye(d*b),ones(1,eta))*Ao5*Hs';
pf5 = kron(eye(d*b),ones(1,eta))*p;

R5 = diag(x5);
S5 = diag(y5);
T5 = diag(z5);

M5 = [zeros(n) -T5         S5;
     T5          zeros(n) -R5;
     -S5         R5          zeros(n)];
 
Am5 = kron(eye(d*b),ones(1,eta))*M5*Ao5*Hs';
pm5 = kron(eye(d*b),ones(1,eta))*M5*p5;

Ab5 = [Af5;Am5];
pb5 = [pf5;pm5];

rankAb5 = rank(Ab5);
dimAb5 = size(Ab5);

fprintf('Ab5 is: \n');
disp(Ab5);
fprintf('pb5 is: \n');
disp(pb5);
fprintf(['\nAb5 has dimensions ' num2str(dimAb5(1)) 'x' num2str(dimAb5(2))...
      ' and rank ' num2str(rankAb5) '.\n']);
fprintf('\nRREF of Ab5 is \n');
disp(rref(Ab5));

% As expected, the perturbations don't really affect the rank because
% forces/moments acting on one vertebra affect the other.

%% Case 6: 3-Vertebra Setup, Extension of Drew's Reference

% Coupling is easy to observe in the first case, but now we want to extend
% our analysis to more vertebrae to see if we can make any generalized
% statements. We add noise, too.

s = 16;
r = 12;
n = 15;
b = 3;
eta = n/b;
d = 3;

Cs = [0 1 0 0 0 0 -1 0 0 0 0 0 0 0 0; % Verticals
      0 0 1 0 0 0 0 -1 0 0 0 0 0 0 0;
      0 0 0 1 0 0 0 0 -1 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 -1 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 -1 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0 -1 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 1 0 0 0 0 -1;
      
      0 0 0 1 0 0 -1 0 0 0 0 0 0 0 0; % Saddles
      0 0 0 1 0 0 0 -1 0 0 0 0 0 0 0;
      0 0 0 0 1 0 -1 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 -1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 -1 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0 -1 0 0;
      0 0 0 0 0 0 0 0 0 1 0 -1 0 0 0;
      0 0 0 0 0 0 0 0 0 1 0 0 -1 0 0];
  
Cr = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0; % Body 1 bars
      1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0;
      1 0 0 -1 0 0 0 0 0 0 0 0 0 0 0;
      1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0;
      
      0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0; % Body 2 bars
      0 0 0 0 0 1 0 -1 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 -1 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 0 -1 0 0 0 0 0;
      
      0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0; % Body 3 bars
      0 0 0 0 0 0 0 0 0 0 1 0 -1 0 0;
      0 0 0 0 0 0 0 0 0 0 1 0 0 -1 0;
      0 0 0 0 0 0 0 0 0 0 1 0 0 0 -1];
C = [Cs;Cr];
Hs = [eye(s) zeros(s,r)];

pinned = zeros(n, 1);
pinned(1:3) = 1;
m_b = 0.8;
m6 = m_b/eta * ones(n, 1);

% Noise vectors by body
xn6p2 = normrnd(0,bar_endpoint/100);
yn6p2 = normrnd(0,bar_endpoint/100);
zn6p2 = normrnd(0,bar_endpoint/100);

xn6a2 = normrnd(0,pi/180);
yn6a2 = normrnd(0,pi/180);
zn6a2 = normrnd(0,pi/180);

xn6p3 = normrnd(0,bar_endpoint/100);
yn6p3 = normrnd(0,bar_endpoint/100);
zn6p3 = normrnd(0,bar_endpoint/100);

xn6a3 = normrnd(0,pi/180);
yn6a3 = normrnd(0,pi/180);
zn6a3 = normrnd(0,pi/180);

xi6 = zeros(b * 6, 1);
xi6(4:6) = [0; pi/2; 0];
xi6(7:12) = [2*bar_endpoint + xn6p2;
            0 + yn6p2;
            0 + zn6p2;
            0 + xn6a2;
            pi/2 + yn6a2;
            0 + zn6a2];
xi6(13:18) = [4*bar_endpoint + xn6p3;
            0 + yn6p3;
            0 + zn6p3;
            0 + xn6a3;
            pi/2 + yn6a3;
            0 + zn6a3];
coordinates6 = get_node_coordinates_3d(a, xi6, debugging);
x6 = coordinates6(1,:)';
y6 = coordinates6(2,:)';
z6 = coordinates6(3,:)';
[px6, py6, pz6] = get_reaction_forces_3d(coordinates6, pinned, m6, g, debugging);
px6 = -px6; py6 = -py6; pz6 = -pz6;
for i=1:n
    pz6(i) = pz6(i) - m6(i)*g;
end
p6 = [px6;py6;pz6];

Ao6 = [C'*diag(C*x6);
       C'*diag(C*y6);
       C'*diag(C*z6)]; % intermediate matrix
Af6 = kron(eye(d*b),ones(1,eta))*Ao6*Hs';
pf6 = kron(eye(d*b),ones(1,eta))*p6;

R6 = diag(x6);
S6 = diag(y6);
T6 = diag(z6);

M6 = [zeros(n) -T6         S6;
     T6          zeros(n) -R6;
     -S6         R6          zeros(n)];
 
Am6 = kron(eye(d*b),ones(1,eta))*M6*Ao6*Hs';
pm6 = kron(eye(d*b),ones(1,eta))*M6*p6;

Ab6 = [Af6;Am6];
pb6 = [pf6;pm6];

rankAb6 = rank(Ab6);
dimAb6 = size(Ab6);

fprintf('Ab6 is: \n');
disp(Ab6);
fprintf('pb6 is: \n');
disp(pb6);
fprintf(['\nAb6 has dimensions ' num2str(dimAb6(1)) 'x' num2str(dimAb6(2))...
      ' and rank ' num2str(rankAb6) '.\n']);
fprintf('\nRREF of Ab6 is \n');
disp(rref(Ab6));

% Comparing to A6, which is the nodal static equilibrium matrix
A6 = [C'*diag(C*x6);
      C'*diag(C*y6);
      C'*diag(C*z6)];
rankA6 = rank(A6);
dimA6 = size(A6);

fprintf('A6 is: \n');
disp(A6);
fprintf('p6 is: \n');
disp(p6);
fprintf(['\nA6 has dimensions ' num2str(dimA6(1)) 'x' num2str(dimA6(2))...
      ' and rank ' num2str(rankA6) '.\n']);
fprintf('\nRREF of A6 is \n');
disp(rref(A6));

% Now this is very interesting... the rank deficiency of the 3 vertebra
% case is now 4, whereas for the 2 vertebra case it was 2. 2/4 of the
% saddle cables between each pair of bodies become force/moment redundant.

% CONJECTURE: for this class of spinal tensegrities, the rigid body
% reformulation will ALWAYS admit solutions for our force density
% optimization problem where the nodal method NEVER will. We can see that
% the nodal matrix retains full rank, just like in the 2-vertebra case.

% CONJECTURE: The rank deficiency of a spinal tensegrity of this form will
% always be 2*(b-1).

% By generality, I expect that adding more vertebrae to the system doesn't
% change the rank relations here, and it seems like the only way we can
% drop rank is by some physically impossible request. These initial results
% suggest that this could be a very promising method.



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

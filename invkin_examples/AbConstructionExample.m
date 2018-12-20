%% Current TODO
% (1) Formally verify moment equation (I'm fairly sure this is correct)
% (2) Prove/disprove full rank of A for this class of tensegrity
% (3) Prove/disprove rank deficiency of Ab for this class of tensegrity
% (4) Write plotting code for 3D systems so we can see what the perturbed
%     states actually look like

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
% system. While this increases s+r, this is not useful for dynamic robot
% systems because they do not add degrees of actuation. Thus, the problem
% will extend into dynamic analysis infinitely for any spinal robots using
% this vertebra structure. Otherwise, we could just trivially transform the
% system into a solvable one by arbitrarily picking points in the middle of
% bars and designating them as nodes.

% The goal of this analysis is to determine the conditions in which we can
% design robots with different formulations of equilibrium matrices such
% that we can abuse rank deficiency and thus produce optimal solutions.

% Stuff to mention: only useful for non-pure tensegrities, systems that
% have frames, etc.

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
C = get_tetrahedral_spine_C_3d(b);
[Hs,~] = get_H(s,r);

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

Ao = get_A(C,x,y,z);
Af = kron(eye(d*b),ones(1,eta))*Ao*Hs';

M = get_M(n,x,y,z);

Ab = get_Ab(d,b,n,Ao,M,Hs);

p = [px;py;pz];
pf = kron(eye(d*b),ones(1,eta))*p;
pm = kron(eye(d*b),ones(1,eta))*M*p;
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

Ao2 = get_A(C,x2,y2,z2);
Af2 = kron(eye(d*b),ones(1,eta))*Ao2*Hs';

M2 = get_M(n,x2,y2,z2);
 
Ab2 = get_Ab(d,b,n,Ao2,M2,Hs);

p2 = [px2;py2;pz2];
pf2 = kron(eye(d*b),ones(1,eta))*p2;
pm2 = kron(eye(d*b),ones(1,eta))*M2*p2;
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

Ao3 = get_A(C,x3,y3,z3);
Af3 = kron(eye(d*b),ones(1,eta))*Ao3*Hs';

M3 = get_M(n,x3,y3,z3);
 
Ab3 = get_Ab(d,b,n,Ao3,M3,Hs);

p3 = [px3;py3;pz3];
pf3 = kron(eye(d*b),ones(1,eta))*p3;
pm3 = kron(eye(d*b),ones(1,eta))*M3*p3;
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

Ao4 = get_A(C,x4,y4,z4);
Af4 = kron(eye(d*b),ones(1,eta))*Ao4*Hs';

M4 = get_M(n,x4,y4,z4);
 
Ab4 = get_Ab(d,b,n,Ao4,M4,Hs);

p4 = [px4;py4;pz4];
pf4 = kron(eye(d*b),ones(1,eta))*p4;
pm4 = kron(eye(d*b),ones(1,eta))*M4*p4;
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

coordinates5 = noisy_coordinate_generator(a,b,bar_endpoint,bar_endpoint/100,pi/180,debugging);

x5 = coordinates5(1,:)';
y5 = coordinates5(2,:)';
z5 = coordinates5(3,:)';
[px5, py5, pz5] = get_reaction_forces_3d(coordinates5, pinned, m, g, debugging);
px5 = -px5; py5 = -py5; pz5 = -pz5;
for i=1:n
    pz5(i) = pz5(i) - m(i)*g;
end
Ao5 = get_A(C,x5,y5,z5);
Af5 = kron(eye(d*b),ones(1,eta))*Ao5*Hs';

M5 = get_M(n,x5,y5,z5);
 
Ab5 = get_Ab(d,b,n,Ao5,M5,Hs);

p5 = [px5;py5;pz5];
pf5 = kron(eye(d*b),ones(1,eta))*p5;
pm5 = kron(eye(d*b),ones(1,eta))*M5*p5;
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

C = get_tetrahedral_spine_C_3d(b);
[Hs,~] = get_H(s,r);

pinned = zeros(n, 1);
pinned(1:3) = 1;
m_b = 0.8;
m6 = m_b/eta * ones(n, 1);

coordinates6 = noisy_coordinate_generator(a,b,bar_endpoint,bar_endpoint/100,pi/180,debugging);

x6 = coordinates6(1,:)';
y6 = coordinates6(2,:)';
z6 = coordinates6(3,:)';
[px6, py6, pz6] = get_reaction_forces_3d(coordinates6, pinned, m6, g, debugging);
px6 = -px6; py6 = -py6; pz6 = -pz6;
for i=1:n
    pz6(i) = pz6(i) - m6(i)*g;
end

Ao6 = get_A(C,x6,y6,z6);
Af6 = kron(eye(d*b),ones(1,eta))*Ao6*Hs';

M6 = get_M(n,x6,y6,z6);
 
Ab6 = get_Ab(d,b,n,Ao6,M6,Hs);

p6 = [px6;py6;pz6];
pf6 = kron(eye(d*b),ones(1,eta))*p6;
pm6 = kron(eye(d*b),ones(1,eta))*M6*p6;
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

%% Case 7: 4-Vertebra Setup
% Just confirming that rank conjectures hold for one more vertebra
s = 24;
r = 16;
n = 20;
b = 4;
eta = n/b;
d = 3;

C = get_tetrahedral_spine_C_3d(b);
[Hs,~] = get_H(s,r);

pinned = zeros(n, 1);
pinned(1:3) = 1;
m_b = 0.8;
m7 = m_b/eta * ones(n, 1);

coordinates7 = noisy_coordinate_generator(a,b,bar_endpoint,bar_endpoint/100,pi/180,debugging);

x7 = coordinates7(1,:)';
y7 = coordinates7(2,:)';
z7 = coordinates7(3,:)';
[px7, py7, pz7] = get_reaction_forces_3d(coordinates7, pinned, m7, g, debugging);
px7 = -px7; py7 = -py7; pz7 = -pz7;
for i=1:n
    pz7(i) = pz7(i) - m7(i)*g;
end

Ao7 = get_A(C,x7,y7,z7);
M7 = get_M(n,x7,y7,z7);

Ab7 = get_Ab(d,b,n,Ao7,M7,Hs);

p7 = [px7;py7;pz7];
pf7 = kron(eye(d*b),ones(1,eta))*p7;
pm7 = kron(eye(d*b),ones(1,eta))*M7*p7;
pb7 = [pf7;pm7];

rankAb7 = rank(Ab7);
dimAb7 = size(Ab7);

fprintf('Ab7 is: \n');
disp(Ab7);
fprintf('pb7 is: \n');
disp(pb7);
fprintf(['\nAb7 has dimensions ' num2str(dimAb7(1)) 'x' num2str(dimAb7(2))...
      ' and rank ' num2str(rankAb7) '.\n']);
fprintf('\nRREF of Ab7 is \n');
disp(rref(Ab7));

A7 = [C'*diag(C*x7);
      C'*diag(C*y7);
      C'*diag(C*z7)];
rankA7 = rank(A7);
dimA7 = size(A7);

fprintf('A7 is: \n');
disp(A7);
fprintf('p7 is: \n');
disp(p7);
fprintf(['\nA7 has dimensions ' num2str(dimA7(1)) 'x' num2str(dimA7(2))...
      ' and rank ' num2str(rankA7) '.\n']);
fprintf('\nRREF of A7 is \n');
disp(rref(A7));

% Confirmed. See case 6 for conjectures on this structure. Also, after 4
% vertebrae, this method will always produce wide matrices by virtue of
% collapsing the problem into body-space, so there are some classes of
% tensegrities that are conclusively solved by this formulation that cannot
% be solved by the nodal one.

%% Extended Analysis
% I now generalize and see if this conjecture holds for much longer spines
% with noise injection in their vertebra positions

vert_max = 50; % Highest vertebra count - these matrices will get massive
bar_endpoint = 0.5; % m
a = [   0,              0,              0;
        bar_endpoint,     0,              -bar_endpoint;
        -bar_endpoint,    0,              -bar_endpoint;
        0,              bar_endpoint,     bar_endpoint;
        0,              -bar_endpoint,    bar_endpoint]';
g = 9.81;

% We care about the following:
% (1) Dimensions of Ab
% (2) Rank of Ab
% (3) Dimensions of A
% (4) Rank of A
for b = 2:vert_max
    % Parameters
    s = 10*(b-1); % SHOULD BE 8, 10 FOR TESTING
    r = 4*b;
    n = 5*b;
    eta = n/b;
    d = 3;
    
    C = get_tetrahedral_spine_C_3d(b);
    [Hs,~] = get_H(s,r);

    pinned = zeros(n, 1);
    pinned(1:3) = 1;
    m_b = 0.8;
    m = m_b/eta * ones(n, 1);

    coordinates = noisy_coordinate_generator(a,b,bar_endpoint,bar_endpoint/100,pi/180,debugging);

    x = coordinates(1,:)';
    y = coordinates(2,:)';
    z = coordinates(3,:)';
    
    A = get_A(C,x,y,z);
    M = get_M(n,x,y,z);
    Ab = get_Ab(d,b,n,A,M,Hs);
    
    rankA = rank(A);
    dimA = size(A);
    rankAb = rank(Ab);
    dimAb = size(Ab);
    
    fprintf(['b = ' num2str(b) ' | dim(A) = ' num2str(dimA) ' | rank(A) = '...
             num2str(rankA) ' | dim(Ab) = ' num2str(dimAb) ' | rank(Ab) = '...
             num2str(rankAb) '\n']);
end

%% Extra: Connectivity Analysis and Graph Coloring
% This is currently commented out
%{

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

% For higher vertebrae systems, this doesn't necessarily immediatley work.
% This is because we actually need to carefully choose a correct basis for
% the null space of the laplacian to pick off node indices. I'm going to
% save this observation for later, but I think we can make something good
% out of it.

%}
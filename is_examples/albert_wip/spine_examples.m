%% Spine Parameters

% --- (1) Tetra Spine ---
tetra = struct;

% source matrix
tetra.Cso = [0 1 0 0 0;
             0 0 1 0 0;
             0 0 0 1 0;
             0 0 0 0 1;
             0 0 0 1 0;
             0 0 0 1 0;
             0 0 0 0 1;
             0 0 0 0 1];

% sink matrix
tetra.Csi = [0 -1 0 0 0;
             0 0 -1 0 0;
             0 0 0 -1 0;
             0 0 0 0 -1;
             0 -1 0 0 0;
             0 0 -1 0 0;
             0 -1 0 0 0;
             0 0 -1 0 0];

% vertebra matrix
tetra.Crv = [1 -1 0 0 0;
             1 0 -1 0 0;
             1 0 0 -1 0;
             1 0 0 0 -1];

tetra.eta = size(tetra.Cso,2); % nodes/vertebra
tetra.sigma = size(tetra.Cso,1); % cables/pair
tetra.rho = size(tetra.Crv,1); % bars/vertebra

%% User-Defined Parameters

% my choice of axis mappings

%    y
%    |                |
%    |                v  gravity
%    |
%    . ------ x
%   /
%  /
% z

% [EDIT] spine parameters
spine = tetra; % spine structure
d = 3; % dimensions
b = 3; % bodies

% [DO NOT EDIT] calculated parameters
[Cso,Csi,Crv,eta,sigma,rho,n,s,r,H] = getParams(spine,b);

% [EDIT] vertebral tranformation matrix. 6 x b.
% [NOTE] first three rows are translations, last three are rotations (rad).
% I choose to define positions absolutely
Uabs = [0  1.5     3;
        0  0.1   -0.1;
        0  0.1   -0.1;
        0 pi/24 -pi/24;
        0 pi/24 -pi/24
        0 pi/24 -pi/24];
U = Uabs;
    
% [EDIT] vertebra mapping matrix. 3 x eta.
V = [0 -sqrt(2) -sqrt(2)  sqrt(2) sqrt(2);
     0 -sqrt(2)  sqrt(2)    0       0;
     0    0        0     -sqrt(2) sqrt(2)];

% [CHOOSE] retrieving coordinates using relative/absolute transformations
% [x,y,z] = getCoordinatesRel(n,U,V);
[x,y,z] = getCoordinatesAbs(n,U,V);

% [EDIT] anchor removal matrix. n x n diag.
% [NOTE] I choose to anchor the first vertebra.
W = kron(eye(d),diag([zeros(1,eta) ones(1,eta*(b-1))]));

% [EDIT] external forces
px = zeros(n,1);
py = -ones(n,1); % placeholder for gravity on nodes
pz = zeros(n,1);

p = W*[px;py;pz]; % [DO NOT EDIT]

% [EDIT] cost functions for optimization
% [NOTE] I choose quadprog and minimum cable tension constraints
minT = 1; % minimum enforced tension

Qn = 2*eye(s+r);
fn = zeros(s+r,1);
Aineqn = -H;
bineq = -ones(s,1) * minT;
opts = optimoptions('quadprog','MaxIterations',100000);
ncostfunc = @(An) quadprog(Qn,fn,Aineqn,bineq,An,p,[],[],[],opts);

Aineqrb = -eye(s);
frb = zeros(s,1);
Qrb = eye(s);
rbcostfunc = @(Ab,beq) quadprog(Qrb,frb,Aineqrb,bineq,Ab,beq,[],[],[],opts);

%% Input Construction

% mandatory inputs
mandInputs = struct(          ...
    'spine',      spine,      ...
    'd',          d,          ...
    'b',          b,          ...
    'eta',        eta,        ...
    'sigma',      sigma,      ...
    'rho',        rho,        ...
    'n',          n,          ...
    's',          s,          ...
    'r',          r,          ...
    'p',          p,          ...
    'x',          x,          ...
    'y',          y,          ...
    'z',          z,          ...
    'U',          U,          ...
    'V',          V,          ...
    'W',          W,          ...
    'H',          H,          ...
    'ncostfunc',  ncostfunc,  ...
    'rbcostfunc', rbcostfunc  ...
);

% optional inputs
optInputs = struct(...
...
);

%% Nodal vs Rigid Body Comparison

% Nodal
q = nodalIK(mandInputs,optInputs);

% Rigid Body
qs = rbIK(mandInputs,optInputs);

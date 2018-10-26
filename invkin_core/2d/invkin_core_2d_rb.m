%% invkin_core_2d_rb.m
% Copyright Andrew P. Sabelhaus 2018

function [f_opt, q_opt, Ab, pb] = invkin_core_2d_rb(x, z, px, pz, C, COMs, s, b, q_min, debugging)
% invkin_core_2d_rb performs a single inverse kinematics computation for a
% 2-dimensional tensegrity structure (robot), using the reformulated
% problem treating the rigid bodies with a force/moment balance and not as
% a node-graph.
%
%   *IMPORANT:* This script currently only is valid under the following
%   conditions:
%
%   d = 2 (this is a 2d script anyway!)
%   b = 2
%   Each rigid body has the same structure (specifically, same number of nodes)
%
% This script follows the force density method for tension networks, which
% calculates a condition for static equilibrium given a specific
% structure, and calculates cable forces / force densities by minimizing
% the total force density in the structure.
%
% (see Schek, Tran, Friesen, etc. for formulation, and my own papers for
% more discussion on 'inverse kinematics'.)
%
% This function REQUIRES that the p vector already be pre-populated with
% reaction forces and gravitational forces, if appropriate.

% Inputs:
%   x, z = the x, and z positions of each node in the structure. Must
%   be the same dimension.
%
%   C = configuration matrix for the structure. See literature for
%   definition.
%
%   px, pz = vectors of external forces, applied to each node (has the same
%   dimension as x, z.)
%
%   COMs = locations of the center of mass of each body. This needs to be
%   pre-computed by the calling function, or else we'd need knowledge of
%   the mass of each rigid body within this script (want to avoid that.)
%   This matrix should be column vectors of each COM, so \in R^{d x b}
%
%   s = number of cables (tension-only members in the system.) This is
%   required to pose the optimization program such that no cables "push."
%
%   b = number of rigid bodies. With fancier algorithms, we could actually
%   do a search through the C matrix (graph search!) to find how many
%   independent cycles there are in the "r" block of it (rods only), but
%   for now, it's easier to just have the caller specify. (Robot designers
%   will know this intuitively: e.g., we've made a robot with 2 vertebrae.)
%
%   q_min = minimum force density for an individual compression member. 
%
%   debugging = set to 1 if you want more information output on the command
%   line.

% Outputs:
%   f_opt = optimal forces (inv kin soln) for ALL members. First s are
%   cables.
%
%   q_opt = optimal force densities, corresponding to f_opt
%
%   Ab, pb = calculated static equilibrium constraint, 
%           provided for debugging if desired.

%% Check if inputs are valid (conformal dimensions of vectors and matrices)

% do this later...

%% First, formulate the constants / parameters for the optimization:

% As aside: we can get the number of nodes and number of members by
n = size(x,1);
r = size(C,1) - s; % rows of C is total number of members, and s is
% provided.

% For later work, generalize: this is the dimensionality of the problem
% (2 dimensional.)
d = 2;

% First, let's create the A matrix.
% We're assuming there are only two rigid bodies.
% That way, we can define the following helper variable that "removes" the
% bar forces from each rigid body.

% Remember, we're looking at something (in the end) like:
% Ab * qs = pb
% where
% Ab = [Af; Am]
% pb = [pf; pm]
% Af = forces from each cable on each rigid body
% pf = sum of external forces in each dimension for each rigid body
%       *NOTE: your calling script should include gravity here already!
% Am = sum of moments from each cable around each rigid body's center of
%       mass
% pm = sum of moments from each set of external reaction forces around each
%       rigid body's center of mass

% Helper matrix:
H_hat = [eye(s), zeros(s, r)];
% ...when multiplied by a vector of all cables, removes the bars.

% Therefore, all the lengths of each cable are:
dx = H_hat * C * x;
dz = H_hat * C * z;

% Since all cables apply forces to both rigid bodies, and assuming tension
% is "negative" here, 
Af = [-dx'; dx'; -dz'; dz'];
% remember, force = length * force density.
  
% Combine the p vector
p = [px; pz];

% Now, we can also calculate pf by "collapsing" it down to per-body instead
% of per-node.
% The number of nodes per rigid body is
eta = n / b;

% a matrix to pattern out the collapsation.
% We could do ones(1, eta) but all Drew's papers use the 'ones' as a column
% vector, so do that here, and transpose it.
collapse = kron(eye(d*b), ones(eta, 1)');

% now, we have pf:
pf = collapse * p;

% Next, we need the moments. This is significantly more complicated.
% The COM for body i is COM(:,i).

% We'll need two vectors for each cable:
% a) from COM of that rigid body to the anchor point for that rigid body,
%       (in paper: v_(h, i)
% b) from the anchor point on this rigid body to its other anchor point
%       (in paper: v_(i, j)

% Though there might be a linear-algebra way to do this next step, a for
% loop is more intuitive.

% Calculate the moments for each cable, for each body.
% Am is rows for each body's x and z sum, \in R^{b x s}
% (since we only have one moment in 2D.)
Am = zeros(b, s);

% For each rigid body,
for g = 1:b
    % for each cable for this body,
    for k = 1:s
        % (note, we know by our assumption of b = 2 that each cable apples
        % a force to both bodies.)
        % From Drew's paper, we have, for body g, and cable k, which
        % connects nodes i and j,
        % (abusing a bit of notation about the COM,)
        % M(g, k) = det( v_{com(g), i}, v_{i, j} ) * q_k
        % 
        % Get each of the two vectors.
        % The anchor points for cable k can be picked out using the C
        % matrix.
        % Remember, the C matrix is graph structured: 1 at the "from" node,
        % and -1 at the "to" node. E.g., anchor1 and anchor2.
        % The == operator checks if each element of the matrix has a
        % certain value.
        from = (C == 1);
        to = (C == -1);
        % HOWEVER, the graph structure of C does NOT account for the
        % "switching" of anchor points between rigid bodies. 
        % For the "first" rigid body, the one with a +1 in C, this makes
        % sense. However, for the "second" rigid body, with the -1 in C,
        % this is flipped: the anchor on this body is actually the
        % *second* node, not the first. 
        % How to check and flip if needed?
        % By our assumption, the "s" rows of C are block structured
        % according to rigid body (e.g., size s x eta each.)
        % So, one easy way is to check if, for the k-th cable,
        % this s-th block has a +1 or -1.
        % That would mean, for an example with b=2, n=8, so eta=4, 
        % "if there is a 1 in to_k, in its columns 5:8, then switch the
        % from and to."
        % A "1" is logical true, and since there can only be at most one
        % anchor point (assuming our C is valid and Drew's math ain't
        % wrong, see the justification for these indices below)
        if( sum( to(k, ((g-1)*eta + 1):(g*eta))) == 1)
            % flip.
            if( debugging )
                disp('Flipping anchors for body')
                g
                disp('and cable')
                k
            end
            temp = from;
            from = to;
            to = temp;
        end
        
        % Remove all but the k-th row of this matrix.
        % By doing so, we've got two vectors that multiply together to get
        % a scalar.
        from_k = from(k, :);
        to_k = to(k,:);

        % We then know the x and z positions of each anchor.
        x_from_k = from_k * x;
        z_from_k = from_k * z;
        x_to_k = to_k * x;
        z_to_k = to_k * z;
        % ...these should all be SCALARS!
        %anchor1 = [C(k, :)*x; C(k, :)*z];
        
        % the COMs matrix: rows are x, z and columns are g.
        % the first vector (COM -> anchor1) is 
        v1 = [x_from_k - COMs(1,g); 
              z_from_k - COMs(2,g)];
        % the second vector (anchor1 -> anchor2) is
        v2 = [x_to_k - x_from_k; 
              z_to_k - z_from_k];
        % The component in the Am matrix is then the determinant
        % (since the q_k vector will be multiplied by it during the
        % optimization.)
        Am(g, k) = det([v1, v2]);
    end
end

% Now, we also need the pm vector.
% The result is a sum, according to body.
pm = zeros(b, 1);
% This sum must be over EACH NODE!
% This is because we need the moment arm for each force!!!
%   ***TO-DO: IS THIS VALID, SINCE GRAVITY IS INCLUDED IN p?
%             ...in other words, since the moment arm for the sum of
%             gravitational forces is zero, since we're summing around the
%             COM... is the sum over gravitational forces still zero since
%             those forces "cancel out" or do we need to do that ourselves?
%             Pretty sure the "axis rule" for moments is such that these will
%             cancel, but can't be sure until I justify it myself.
for g = 1:b
    % For the g-th rigid body, we will pull out the coordinates for each of
    % its nodes.
    % Nodes go from ((g-1)*eta + 1) : g*eta, 
    % example, 2 bodies 8 nodes, body 1 is 1:4, body 2 is 5:8.
    x_g = x( (g-1)*eta + 1 : g*eta );
    z_g = z( (g-1)*eta + 1 : g*eta );
    % Similarly, pull out the forces from the p vectors (saves us indexing
    % trouble later.)
    px_g = px( (g-1)*eta + 1 : g*eta );
    pz_g = pz( (g-1)*eta + 1 : g*eta );
    % We'll form the total force vector, knowing that the direction of the
    % p_x forces is [1; 0] and the direction of the p_z forces is [0; 1].
    % This has a slightly different shape than "p" above, in order to work
    % well with the moment determinant. Specifically,
    % p_g \in R^{d, eta}
    p_g = [px_g'; pz_g'];
    % The moment due to node m (from m=1:eta) is similar to above, but we
    % don't need to factor out the force itself (can leave in, unlike the
    % force densitites q above, so we're really doing r x F here (in 2D).
    % M(g, m) = det( v_{com(g), m}, p_g(m) )
    for m=1:eta
        % This node's vector from the COM for its body is
        v1 = [x_g(m) - COMs(1,g);
              z_g(m) - COMs(2,g)];
        % This node's contribution uses the p_g vector for node m.
        pm(g) = pm(g) + det([v1, p_g(:,m)]);
    end
end

% finally, can concatenate the force balance with the moment balance.
Ab = [Af; Am];
pb = [pf; pm];

if debugging
    disp('Size of x, r, C, Ab, pb is:');
    n
    r
    size(C)
    size(Ab)
    size(pb)
    disp('Values of Ab and pb are:');
    Ab
    pb
end

% Since we assume that the first s rows of C are for the cables, make the
% constraint for the min force density on those cables:
q_min_vector = q_min * [ones(s,1)];

%% Solve the optimization problem

% for now, we're going to let quadprog do the relaxation from equality
% constrained to inequality constrained.

if debugging
    disp('Solving the inverse kinematics problem for a 3D tensegrity structure...');
end

% quadprog's arguments are
% 0.5 * H (quadratic term), f (linear term), A_ineq, b_ineq, A_eq, b_eq

% Our problem is
% min qs'qs s.t. Ab*qs=pb, q - mintensionvec \geq 0

% no weighting any members different than others
% the 2 is unnecessary here since we only care about argmax
%H = 2 * eye(s + r);

% Our H is only in s x s now
H = eye(s);

% no linear term
f = [];

% the inequality constraints are in the form of A_eq * q <= b_eq,
% so we need q >= mintensionvec becomes -q <= -mintensionvec
A_ineq = - eye(s);
b_ineq = - q_min_vector;

% equality constraints are A q = p
A_eq = Ab;
b_eq = pb;

if debugging
    disp('Optimization parameters are:');
    H
    f
    A_eq
    b_eq
end

% Call quadprog
[q_opt] = quadprog(H, f, A_ineq, b_ineq, A_eq, b_eq);

if debugging
    disp('Optimal force densities are:');
    q_opt
end

%% Calculate the cable forces from the force densities.

% % Forces are density * length, and lengths for each member in each
% % direction are
% dx = C * x;
% dz = C * z;

% so the lengths of each cable are the euclidean norm of each 2-vector.
% re-organize:
length_vecs = [dx'; dz'];
if debugging
    disp('Member length vectors are:');
    length_vecs
end

% the scalar lengths are then the 2-norm (euclidean) for each column, which
% is
lengths = vecnorm(length_vecs);
if debugging
    disp('Member lengths (scalar) are:');
    lengths
end

% Then, calculate the optimal forces. (results should be same dim as
% q_opt.)
% again, probably some nice linear algebra here (TO-DO: efficiency),
% but looping is intuitive and I'm busy.
f_opt = zeros(size(q_opt));
for k=1:s
    f_opt(k) = lengths(k) * q_opt(k);
end
%f_opt = lengths * q_opt;

if debugging
    disp('Optimal forces on cables are:');
    f_opt
end

end















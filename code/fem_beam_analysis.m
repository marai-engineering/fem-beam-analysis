
% FEM Beam Analysis of a Bent Beam Structure
% Marai Abed Alrahman
% MATLAB implementation using Euler-Bernoulli beam elements
% Includes:
% - static analysis
% - reaction forces
% - deflection, rotation, and bending moment plots
% - midpoint deflections
% - modal analysis (natural frequencies)
clear; clc; close all;

%% --------------------- INPUT DATA ---------------------
% Geometry (m)
a = 0.220;
b = 0.540;
c = 0.730;          % total length

% Cross-sections
d1 = 0.031;         % [m]
d2 = 2*d1;          % [m]

% Material
E2  = 30e9;         % [Pa]
E1  = 4*E2;         % [Pa]  % small-diameter part = 4E2
rho = 2500;         % [kg/m^3]

% Loads
p1 = 3000;          % [N/m]  (uniform, upward, on elems 1 and 2)
F1 = 2000;          % [N]    (concentrated vertical force at node 3)
M1 = 900;           % [Nm]   (concentrated torque at node 2)

%% --------------------- MESH & SECTIONS ----------------
% Nodes along the beam: 1 -- 2 -- 3 -- 4
x     = [0; a; b; c];   % global x-positions of nodes
nNode = numel(x);
nDof  = 2*nNode;        % [w1,phi1,w2,phi2,w3,phi3,w4,phi4]

% Elements: 3 beam elements between the 4 nodes
elem  = [1 2;
         2 3;
         3 4];
nElem = size(elem,1);

% Element lengths (matching theory: L1=a, L2=b-a, L3=c-b)
L1 = a;
L2 = b - a;
L3 = c - b;
L  = [L1; L2; L3];

% Second moments of area and EI
I1 = pi*d1^4/64;
I2 = pi*d2^4/64;
EI = [E1*I1;   E2*I2;   E2*I2];

% Areas for mass matrix
A1 = pi*d1^2/4;
A2 = pi*d2^2/4;
A  = [A1; A2; A2];

%% --------------------- GLOBAL MATRICES ----------------
K = zeros(nDof);         % global stiffness
M = zeros(nDof);         % global mass
F = zeros(nDof,1);       % global load vector

for e = 1:nElem
    Le  = L(e);
    EIe = EI(e);
    Ae  = A(e);

    i = elem(e,1);
    j = elem(e,2);
    edofs = [2*i-1 2*i 2*j-1 2*j];   % [wi,phii,wj,phij]

    % Beam stiffness matrix (Euler–Bernoulli, 2DOF/node)
    ke = (EIe/Le^3) * [ ...
        12      6*Le   -12      6*Le;
        6*Le  4*Le^2   -6*Le  2*Le^2;
       -12     -6*Le    12     -6*Le;
        6*Le  2*Le^2   -6*Le  4*Le^2 ];

    % Consistent mass matrix
    me = (rho*Ae*Le/420) * [ ...
        156     22*Le    54    -13*Le;
        22*Le  4*Le^2  13*Le   -3*Le^2;
        54     13*Le   156    -22*Le;
       -13*Le -3*Le^2 -22*Le   4*Le^2];

    % Equivalent nodal loads from p1 on elements 1 and 2
    feq = zeros(4,1);
    if e <= 2
        % standard equivalent nodal loads for a uniform q = p1
        feq = p1*Le/2 * [1; Le/6; 1; -Le/6];
    end

    % Assembly
    K(edofs,edofs) = K(edofs,edofs) + ke;
    M(edofs,edofs) = M(edofs,edofs) + me;
    F(edofs)       = F(edofs)       + feq;
end

%% --------------------- CONCENTRATED LOADS ----------------
% vect(0,0,0,M1,F1,0,0,0) -> DOF4 = M1, DOF5 = F1
F(4) = F(4) + M1;    % torque at node 2 (phi2)
F(5) = F(5) + F1;    % vertical force at node 3 (w3)

%% --------------------- BOUNDARY CONDITIONS ----------------
% Fixed DOFs are 3, 7, 8:
% DOF map: [1 2 3 4 5 6 7 8] = [w1,phi1,w2,phi2,w3,phi3,w4,phi4]
fixed   = [3 7 8];
allDofs = 1:nDof;
free    = setdiff(allDofs,fixed);

%% --------------------- 1–2. STATIC SOLUTION --------------------
U       = zeros(nDof,1);
U(free) = K(free,free) \ F(free);   % Kff*Uf = Ff
R       = K*U - F;                  % reactions

fprintf('\n=== STATIC RESULTS ===\n');
fprintf('Nodal displacements (w in mm, phi in rad):\n');
for i = 1:nNode
    w   = 1e3*U(2*i-1);      % [mm]
    phi = U(2*i);            % [rad]
    fprintf('Node %d: w = %+8.4f mm   phi = %+9.4e rad\n', i, w, phi);
end

fprintf('\nReaction forces / moments at fixed DOFs:\n');
for k = fixed
    node = ceil(k/2);
    if mod(k,2)==1          % odd DOF -> w
        fprintf('DOF %d (node %d, w ): R = %+10.3f N\n',  k, node, R(k));
    else                    % even DOF -> phi
        fprintf('DOF %d (node %d,phi): R = %+10.3f Nm\n', k, node, R(k));
    end
end

%% --------------------- 3. FIELD PLOTS (w, phi, M_b) -------------
nPlot = 80;
X = []; W = []; PHI = []; Mdiag = [];

for e = 1:nElem
    i   = elem(e,1);
    j   = elem(e,2);
    Le  = L(e);
    EIe = EI(e);

    w1  = U(2*i-1);   phi1 = U(2*i);
    w2  = U(2*j-1);   phi2 = U(2*j);

    xi  = linspace(0,1,nPlot);      % local 0..1

    % Hermite shape functions
    N1  = 1 - 3*xi.^2 + 2*xi.^3;
    N2  = Le*(xi - 2*xi.^2 + xi.^3);
    N3  = 3*xi.^2 - 2*xi.^3;
    N4  = Le*(-xi.^2 + xi.^3);

    wloc = N1.*w1 + N2.*phi1 + N3.*w2 + N4.*phi2;

    % Derivatives wrt xi (for rotation)
    N1p = -6*xi + 6*xi.^2;
    N2p = Le*(1 - 4*xi + 3*xi.^2);
    N3p =  6*xi - 6*xi.^2;
    N4p = Le*(-2*xi + 3*xi.^2);

    % rotation = dw/dx = (1/Le)*dw/dxi
    philoc = (1/Le)*(N1p.*w1 + N2p.*phi1 + N3p.*w2 + N4p.*phi2);

    % curvature and bending moment (sign chosen so M positive in sagging)
    d2 = (6/Le^2)*(w1*(2*xi-1) + w2*(1-2*xi)) ...
       + (2/Le)*(phi1*(3*xi-2) + phi2*(3*xi-1));
    Mloc = -EIe * d2;

    Xe = x(i) + xi*(x(j)-x(i));  % map to global x

    X     = [X Xe];
    W     = [W wloc];
    PHI   = [PHI philoc];
    Mdiag = [Mdiag Mloc];
end

figure;
subplot(3,1,1)
plot(X,1e3*W,'LineWidth',2); grid on;
ylabel('w(x) [mm]');
title('Deflection');

subplot(3,1,2)
plot(X,PHI,'LineWidth',2); grid on;
ylabel('\phi_z(x) [rad]');
title('Rotation');

subplot(3,1,3)
plot(X,Mdiag/1e3,'LineWidth',2); grid on;
ylabel('M_b(x) [kNm]');
xlabel('x [m]');
title('Bending moment');

%% --------------------- 4. DEFLECTION AT MID–ELEMENT POINTS ----
% midpoints at x = a/2, (a+b)/2, (b+c)/2 (one per element)
w_mid = zeros(3,1);

for e = 1:nElem
    i   = elem(e,1);
    j   = elem(e,2);
    Le  = L(e);

    w1  = U(2*i-1);   phi1 = U(2*i);
    w2  = U(2*j-1);   phi2 = U(2*j);

    xi  = 0.5;        % element midpoint

    N1  = 1 - 3*xi^2 + 2*xi^3;
    N2  = Le*(xi - 2*xi^2 + xi^3);
    N3  = 3*xi^2 - 2*xi^3;
    N4  = Le*(-xi^2 + xi^3);

    w_mid(e) = N1*w1 + N2*phi1 + N3*w2 + N4*phi2;
end

fprintf('\nDeflection at midpoints (for table 3):\n');
fprintf('x = a/2       -> w = %+8.4f mm\n', 1e3*w_mid(1));
fprintf('x = (a+b)/2   -> w = %+8.4f mm\n', 1e3*w_mid(2));
fprintf('x = (b+c)/2   -> w = %+8.4f mm\n', 1e3*w_mid(3));

%% --------------------- 5. NATURAL FREQUENCIES -----------------
% Free–vibration eigenproblem: Kff*phi = w^2 * Mff * phi (unloaded)
Kff = K(free,free);
Mff = M(free,free);

[Phi_modes,Omega2] = eig(Kff,Mff);           % generalized eigenproblem
omega = sqrt(real(diag(Omega2)));            % [rad/s]
[omega,idx] = sort(omega);                   % sort
freq = omega/(2*pi);                         % [Hz]

nModes = min(3,numel(freq));
fprintf('\n=== FIRST %d NATURAL FREQUENCIES (UNLOADED) ===\n', nModes);
for m = 1:nModes
    fprintf('Mode %d:  omega = %8.3f rad/s   f = %8.3f Hz\n', ...
        m, omega(m), freq(m));
end

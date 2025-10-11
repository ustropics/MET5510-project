global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These can be found in Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% beginning on slides 5 - 6

%% Initial variables for complex and wave number
cplx = sqrt(-1); % imaginary unit for complex number operations
m0 = 7; % wave number (set to 7 as in examples)

%% grid point parameterization 
jj = 50; % number of latitude grid points
kk = 50; % number of height grid points
ll = (jj-1)*(kk+1); % total number of linear indices for interior points
Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)

%% Coriolis, beta, and gravity calculations/constants
f0 = 2 * (2*pi/86400) * sin(pi/4); % coriolis parameter at 45° latitude (s^-1)
beta = 2 * (2 * pi/86400) * cos(pi/4)/(6.37*1.e6); % beta parameter (meridional PV gradient, m^-1 s^-1) - non-zero for beta-plane
gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures, height, and additional domain gridding
Theta0 = 300; % reference potential temperature (K)
delta_Theta0 = 30; % vertical potential temperature difference (K) - set to 30 as per model
HH = 10000; % scale height (m)

Ly = 6.37 * 1.0e6*50 * pi/180; % meridional domain length (50° latitude range, m)
dy = Ly/jj; % meridional grid spacing (in meters)
dz = HH/kk; % vertical grid spacing (in meters)

%% Set brunt-vaisala frequency
NN2 = gg * delta_Theta0/(HH*Theta0); % brunt-vaisala frequency squared (s^-2)

%% Parameters for Modified Hoskins-West Model
DeltaT_hor = 60; % horizontal temperature difference (K)
mu = 0.5; % mu parameter (0.5 or 1)
U0 = 0; % mean zonal wind offset (set to 0)
L_d = sqrt(NN2) * HH / f0; % Rossby deformation radius
m_y = 2 * pi / Ly; % meridional wavenumber
gamma = m_y * L_d; % gamma for sinh term
y0 = Ly / 2; % center y
Lambda = gg * DeltaT_hor / (f0 * Theta0 * Ly); % shear rate

%% Set Ubar with Modified Hoskins-West formula
Ubar = zeros(jj+1,kk+1); % initialize mean zonal wind field
for j = 1:jj+1
    y = (j-1)*dy; % y coordinate
    for k = 1:kk+1
        z = (k-1)*dz; % z coordinate
        term_linear = z / HH;
        term_quadratic = -mu / 2 * (z / HH)^2;
        if abs(gamma) > 1e-10  % avoid division by zero
            term_sinh = sinh(gamma * z / HH) / sinh(gamma) * cos(m_y * (y - y0));
        else
            term_sinh = 0;
        end
        Ubar(j,k) = Lambda * HH * (term_linear + term_quadratic + term_sinh) - U0;
    end
end

%% Initialize BPVy = zeros
BPVy = zeros(jj+1,kk+1); % initialize beta-plane PV gradient

%% Set BPVy for interior points using analytical formula
for j = 2:jj
    y = (j-1)*dy;
    for k = 2:kk
        z = (k-1)*dz;
        if abs(gamma) > 1e-10
            sinh_term = sinh(gamma * z / HH) / sinh(gamma);
        else
            sinh_term = 0;
        end
        BPVy(j,k) = beta - 2 * Lambda * HH * m_y^3 * sinh_term * sin(m_y * (y - y0));
    end
end

%% This is just for testing
[j,k] = l2jk(98); % convert linear index 98 to 2D indices (j,k)
y = jk2l(j,k); % convert back to linear index for testing

%% Initialize vectors for PV and its advections
XV = zeros(ll,1); % streamfunction vector
QV = zeros(ll,1); % PV vector
xQVadv = zeros(ll,1); % zonal advection of PV
yQVadv = zeros(ll,1); % meridional advection of PV

%% Initialize matrices
B = zeros(ll,ll); % matrix for PV inversion
C = zeros(ll, ll); % matrix for zonal PV advection
D = zeros(ll,ll); %  matrix for meridional PV advection (beta term = constant)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP & EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main loop to construct matrices B, C, and D
% This is from page 16 from Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% (Four vectors involved in the model)
for l0 = 1:ll
    XV = zeros(ll,1); % reset streamfunction vector
    XV(l0) = 1;  % set 10m to level 1

    % B MATRIX
    QV = stream2pv(XV);
    B(:, l0) = QV(:);

    % C MATRIX
    xQVadv = stream2xPVadv(QV);
    C(:, l0) = xQVadv(:);

    % D MATRIX
    yQVadv = stream2yPVadv(XV);
    D(:,l0) = yQVadv(:);

end

%% Start eigenvector and eigenvalue part
% Linearized PV equation in matrix form for eigenvalue problem
A = B^(-1) * (C+D);

[eigVec, eigValm] = eig(A); % eigVec = eigenvectors
eigVal = diag(eigValm); % eigValm = eigenvalue matrix

% test values in descending order
[test, sortdx] = sort(real(eigVal), 'descend');


% this is a solution to the real part
eigVal2 = eigVal(sortdx);
eigVec2 = eigVec(:, sortdx);

% test values in ascending order
[test, sortdx] = sort(real(cplx*eigVal), 'ascend');

% this is a solution for the complex/imaginary part
eigVal3 = eigVal(sortdx);
eigVec3 = eigVec(:,sortdx);

toc

save('HoskinsWest_wave.mat') % save data file for Modified Hoskins-West Model
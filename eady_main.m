global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These can be found in Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% beginning on slides 5 - 6

%% Initial variables for complex and wave number
cplx = sqrt(-1); % imaginary unit for complex number operations
m0 = 8; % CHANGED: wave number (set to 8 as in PDF page 3 for Eady shortwave; was 2 for Rossby longwave)

%% grid point parameterization 
jj = 50; % number of latitude grid points
kk = 50; % number of height grid points
ll = (jj-1)*(kk+1); % total number of linear indices for interior points
Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)

%% Coriolis, beta, and gravity calculations/constants
f0 = 2 * (2*pi/86400) * sin(pi/4); % coriolis parameter at 45° latitude (s^-1)
beta = 0; % CHANGED: set to 0 for f-plane (Eady model; was computed with cos term for beta-plane/Rossby)
gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures, height, and additional domain gridding
Theta0 = 300; % reference potential temperature (K)
delta_Theta0 = 30; % potential temperature difference across dmain (K) -- this is VERTICAL for N^2
HH = 10000; % scale height (m)

Ly = 6.37 * 1.0e6*50 * pi/180; % meridional domain length (50° latitude range, m)
dy = Ly/jj; % meridional grid spacing (in meters)
dz = HH/kk; % vertical grid spacing (in meters)

%% ADDED: Horizontal temperature difference for meridional gradient (Eady model)
DeltaT_hor = 50; % K (as in PDF page 3; adjust to 30, 40, or 60 as needed)

%% ADDED: Compute shear Lambda from thermal wind
Lambda = gg * DeltaT_hor / (f0 * Theta0 * Ly); % s^-1 (du/dz)

%% Set brunt-vaisala frequency, ubar, and beta-plane (constant)
NN2 = gg * delta_Theta0/(HH*Theta0); % brunt-vaisala frequency squared (s^-2)
Ubar = zeros(jj+1,kk+1); % initialize mean zonal wind field
% CHANGED: Set to linear shear for Eady (u = Lambda * z; z=0 at bottom)
for j = 1:jj+1
    for k = 1:kk+1
        Ubar(j,k) = Lambda * ((k-1)*dz); % m/s (increases with height; ~0 at bottom, ~30-35 at top)
    end
end
BPVy = zeros(jj+1,kk+1); % initialize beta-plane PV gradient

%% This is just for testing
[j,k] = l2jk(98); % convert linear index 98 to 2D indices (j,k)
y = jk2l(j,k); % convert back to linear index for testing

%% set BPVy to beta for interior points
for k = 2 : kk
    for j = 2:jj
        BPVy(j,k) = beta; % assign constant beta to the PV gradient (will be 0 for Eady)
    end
end

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

save('eady_wave.mat') % CHANGED: save to new file for Eady model (was 'Rossby_wave_2.mat')
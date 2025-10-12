%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_main.m

% Description: Main script for setting up and solving the linear 
% quasi-geostrophic (QG) model for Rossby waves, including grid 
% initialization, mean flow (Ubar) and PV gradient (BPVy) setup, matrix 
% construction (B, C, D), and eigenvalue/eigenvector computation for 
% wave stability analysis.

% Functions used: 
% - l2jk: Converts linear index to 2D indices (j, k).
% - jk2l: Converts 2D indices (j, k) to linear index.
% - stream2pv: Computes potential vorticity from streamfunction.
% - stream2xPVadv: Calculates zonal advection of potential vorticity.
% - stream2yPVadv: Computes meridional advection of potential vorticity.
% - eig: Performs eigenvalue decomposition.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('functions')

global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These can be found in Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% beginning on slides 5 - 6

%% Initial variables for complex and wave number
cplx = sqrt(-1); % imaginary unit for complex number operations
m0 = 2; % this is wave number (0-3 is longwave/planetary, 3+ is shortwave)

%% grid point parameterization 
jj = 50; % number of latitude grid points
kk = 50; % number of height grid points
ll = (jj-1)*(kk+1); % total number of linear indices for interior points
Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)

%% Coriolis, beta, and gravity calculations/constants
f0 = 2 * (2*pi/86400) * sin(pi/4); % coriolis parameter at 45° latitude (s^-1)
beta = 2 * (2 * pi/86400) * cos(pi/4)/(6.37*1.e6); % beta parameter (meridional PV gradient, m^-1 s^-1)
gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures, height, and additional domain gridding
Theta0 = 300; % reference potential temperature (K)
delta_Theta0 = 30; % potential temperature difference across dmain (K)
HH = 10000; % scale height (m)

Ly = 6.37 * 1.0e6*50 * pi/180; % meridional domain length (50° latitude range, m)
dy = Ly/jj; % meridional grid spacing (in meters)
dz = HH/kk; % vertical grid spacing (in meters)

%% Set brunt-vaisala frequency, ubar, and beta-plane (constant)
NN2 = gg * delta_Theta0/(HH*Theta0); % brunt-vaisala frequency squared (s^-2)
Ubar = zeros(jj+1,kk+1); % initialize mean zonal wind field (constant 10 m s^-1)
Ubar(:,:) = Ubar(:,:) + 10; % set uniform zonal wind speed
BPVy = zeros(jj+1,kk+1); % initialize beta-plane PV gradient

%% This is just for testing
[j,k] = l2jk(98); % convert linear index 98 to 2D indices (j,k)
y = jk2l(j,k); % convert back to linear index for testing

%% set BPVy to beta for interior points
for k = 2 : kk
    for j = 2:jj
        BPVy(j,k) = beta; % assign constant beta to the PV gradient
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP & EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% create data folder if it doesn't exist
if ~exist('data', 'dir')
    mkdir('data');
end

% create filename and save it
filename = ['data/rossby_wave_', num2str(m0), '.mat'];
save(filename);
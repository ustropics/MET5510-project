%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: hwm_main.m

% Description: Script for setting up the Modified Hoskins-West Model in the 
% quasi-geostrophic framework, initializing grid parameters, computing the 
% mean zonal wind (Ubar) and PV gradient (BPVy) with analytical formulas, 
% constructing matrices (B, C, D) for PV inversion and advection, and 
% solving the eigenvalue problem to determine wave stability, saving results 
% to 'data/HoskinsWest_wave.mat'.

% Functions used: 
% - l2jk: Converts linear index to 2D indices (j, k).
% - jk2l: Converts 2D indices (j, k) to linear index.
% - stream2pv: Computes potential vorticity from streamfunction.
% - stream2xPVadv: Calculates zonal advection of potential vorticity.
% - stream2yPVadv: Computes meridional advection of potential vorticity.
% - eig: Performs eigenvalue decomposition.
% - ubarCalc: Computes mean zonal wind field.
% - BPVyCalc: Computes beta-plane PV gradient.
% - matricesBCD: Constructs matrices B, C, D for PV inversion and advection.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx
tic

addpath('functions'); % add functions folder
addpath('config'); % add config folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These can be found in Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% beginning on slides 5 - 6

%% Load constants from hwm_config.m and assign global variables
params = hwm_config();

% get directories
hwm_data_dir = params.hwm_data_dir;
hwm_data_filename = params.hwm_data_filename;
hwm_data = fullfile(hwm_data_dir, hwm_data_filename);

cplx = params.cplx; % imaginary unit for complex number operations
m0 = params.m0; % wave number (set to 7 as in examples)
jj = params.jj; % number of latitude grid points
kk = params.kk; % number of height grid points
ll = params.ll; % total number of linear indices for interior points
Lx = params.Lx; % zonal domain length at 45° latitude (m)
f0 = params.f0; % coriolis parameter at 45° latitude (s^-1)
beta = params.beta; % beta parameter (meridional PV gradient, m^-1 s^-1)
gg = params.gg; % gravity (m s^-2)
Theta0 = params.Theta0; % reference potential temperature (K)
delta_Theta0 = params.delta_Theta0; % vertical potential temperature difference (K)
HH = params.HH; % scale height (m)
Ly = params.Ly; % meridional domain length (50° latitude range, m)
dy = params.dy; % meridional grid spacing (in meters)
dz = params.dz; % vertical grid spacing (in meters)
NN2 = params.NN2; % brunt-vaisala frequency squared (s^-2)
DeltaT_hor = params.DeltaT_hor; % horizontal temperature difference (K)
mu = params.mu; % mu parameter (0.5 or 1)
U0 = params.U0; % mean zonal wind offset (set to 0)
L_d = params.L_d; % Rossby deformation radius
m_y = params.m_y; % meridional wavenumber
gamma = params.gamma; % gamma for sinh term
y0 = params.y0; % center y
Lambda = params.Lambda; % shear rate

%% Set Ubar with Modified Hoskins-West formula
Ubar = hwm_ubar(params);

%% Initialize BPVy and set interior points using analytical formula
BPVy = hwm_BPVyCalc(params);

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
D = zeros(ll,ll); % matrix for meridional PV advection (beta term = constant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP & EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main loop to construct matrices B, C, and D
[B, C, D] = matricesBCD(ll, @stream2pv, @stream2xPVadv, @stream2yPVadv);

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

%% Create data folder if it doesn't exist and save data file
if ~exist(hwm_data_dir, 'dir')
    mkdir(hwm_data_dir);
end

% save it
save(hwm_data, 'Ubar', 'BPVy', 'eigVec', 'eigVal', 'eigVec2', 'eigVal2', ...
    'eigVec3', 'eigVal3', 'jj', 'kk', 'll', 'NN2', 'f0', 'dy', 'dz', 'm0', ...
    'Lx', 'beta', 'cplx', 'gg', 'Theta0', 'HH', 'Ly');
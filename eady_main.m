%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: eady_main.m

% DESCRIPTION: Main script for the Eady wave model in the quasi-geostrophic
% framework. Initializes parameters, computes Ubar and BPVy, constructs 
% matrices for PV inversion and advection, solves the eigenvalue problem for 
% wave stability, and saves results.

% SCRIPTS:
% - eady_config.m: Loads model parameters

% FUNCTIONS:
% - eady_ubar: Initializes mean zonal wind field
% - bpvy: Initializes potential vorticity gradient field
% - matrices: Constructs matrices B, C, D for PV inversion and advection
% - eigen: Solves eigenvalue problem and sorts results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll HH ZZ YY BPVy NN2 f0 gg dy m0 dz Lx Theta0 dTbar ... 
    Ubar Ubar0 beta cplx model eady_data

addpath('functions'); % add functions folder
addpath('config'); % add config folder
addpath('plots') % add plots folder

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load parameters from eady_config.m
params = eady_config();

% Assign global variables
model = params.model; % model name
cplx = params.cplx; % imaginary unit
m0 = params.m0; % wave number
jj = params.jj; % number of latitude grid points
kk = params.kk; % number of height grid points
ll = params.ll; % total number of linear indices
Lx = params.Lx; % zonal domain length (m)
f0 = params.f0; % Coriolis parameter (s^-1)
beta = params.beta; % beta parameter (m^-1 s^-1)
gg = params.gg; % gravity (m s^-2)
Theta0 = params.Theta0; % reference potential temperature (K)
HH = params.HH; % scale height (m)
dy = params.dy; % meridional grid spacing (m)
dz = params.dz; % vertical grid spacing (m)
NN2 = params.NN2; % Brunt-Vaisala frequency squared (s^-2)
ZZ = params.ZZ; % vertical grid coordinates (m)
YY = params.YY; % meridional grid coordinates (m)
dTbar = params.dTbar; % potential temperature gradient (K)
Ubar0 = params.Ubar0; % base zonal wind speed (m/s)
eady_data = fullfile(params.eady_data_dir, params.eady_data_filename);

%% Initialize Ubar and BPVy
Ubar = eady_ubar(jj, kk, gg, f0, Theta0, dTbar, params.Ly, ZZ, Ubar0);
BPVy = eady_bpvy(jj, kk, beta, Ubar, dy, f0, NN2, dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% MATRIX CONSTRUCTION & EIGENVALUES %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct matrices B, C, D
[B, C, D] = matrices(ll, @stream2pv, @stream2xPVadv, @stream2yPVadv);

%% Solve eigenvalue problem
[eigVec, eigVal, eigVec2, eigVal2, eigVec3, eigVal3] = eigen(B, C, D, cplx);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(params.eady_data_dir, 'dir')
    mkdir(params.eady_data_dir);
end
save(eady_data, 'Ubar', 'BPVy', 'eigVec', 'eigVal', 'eigVec2', 'eigVal2', ...
    'eigVec3', 'eigVal3', 'jj', 'kk', 'll', 'NN2', 'f0', 'dy', 'dz', 'm0', ...
    'Lx', 'beta', 'cplx', 'gg', 'Theta0', 'HH', 'B', 'C', 'D', 'model');

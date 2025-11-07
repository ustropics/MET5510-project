%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: main.m

% DESCRIPTION: Main script for the Hoskins-West modified Eady-type model in the
% quasi-geostrophic framework. Initializes parameters, computes mean zonal wind
% (Ubar) and PV gradient (BPVy), constructs matrices for PV inversion and
% advection, solves the eigenvalue problem for wave stability, and saves results.

% This version uses the Hoskins–West Modified Eady-type mean flow from
% Sample_basic_flows.pdf (p. 9):
%
%   U(y,z) = (g / (f0 * Theta0)) * (H * ΔT / Ly) *
%             [ (z/H) - (μ/2)*(z/H)^2
%               + sinh(2πLy/Lx * z/H)/sinh(2πLy/Lx) * cos(2π(y - y0)/Lx) ]
%
% where μ = 0 → Eady model; μ = 0.5 or 1 → Hoskins–West Eady-type.

% SCRIPTS:
% - cfg.m: Loads model parameters

% FUNCTIONS:
% - ubar: Initializes mean zonal wind field
% - bpvy: Initializes potential vorticity gradient field
% - matrices: Constructs matrices B, C, D for PV inversion and advection
% - eigen: Solves eigenvalue problem and sorts results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll HH ZZ YY BPVy NN2 f0 gg dy m0 dz Lx Theta0 dTbar N mu ... 
    Ubar beta cplx model Lr Ly data

fprintf('Starting to compile data for Quasi-Geostrophic Cyclogenesis, Advection, and Instability (QG-CAI) model...\nthis might take ~1 minute (BRB cig break)...')
tic

addpath('functions'); % add functions folder
addpath('plots') % add plots folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load parameters from config.m
params = cfg();

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
mu = params.mu; % curvature parameter
N = params.N; % Brunt-Vaisala frequency (s^-1)
Lr = params.Lr; % Rossby radius of deformation (m)
Ly = params.Ly; % meridional domain length (50° latitude range, m)
data = fullfile(params.data_dir, params.data_filename);

%% Initialize Ubar and BPVy
Ubar = ubar(jj, kk, gg, f0, Theta0, dTbar, HH, Ly, ZZ, YY, mu, Lr);
BPVy = bpvy(model, jj, kk, beta, Ubar, dy, f0, NN2, dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% MATRIX CONSTRUCTION & EIGENVALUES %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct matrices B, C, D
[B, C, D] = matrices(ll, @stream2pv, @stream2xPVadv, @stream2yPVadv);

%% Solve eigenvalue problem
A = B^(-1) * (C + D);
[eigVec, eigValm] = eig(A); % eigVec = eigenvectors, eigValm = eigenvalue matrix
eigVal = diag(eigValm);

% Sort by descending real part
[~, sortdx] = sort(real(eigVal), 'descend');
eigVal2 = eigVal(sortdx);
eigVec2 = eigVec(:, sortdx);

% Sort by ascending real part of imaginary eigenvalue
[~, sortdx] = sort(real(cplx * eigVal), 'ascend');
eigVal3 = eigVal(sortdx);
eigVec3 = eigVec(:, sortdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(params.data_dir, 'dir')
    mkdir(params.data_dir);
end
save(data, 'Ubar', 'BPVy', 'eigVec', 'eigVal', 'eigVec2', 'eigVal2', ...
    'eigVec3', 'eigVal3', 'jj', 'kk', 'll', 'NN2', 'f0', 'dy', 'dz', 'm0', ...
    'Lx', 'beta', 'cplx', 'gg', 'Theta0', 'HH', 'B', 'C', 'D', 'model');

toc

disp_str = ['All data computed and saved to ', params.data_dir, params.data_filename];
disp(disp_str);
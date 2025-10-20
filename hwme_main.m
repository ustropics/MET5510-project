%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: hwme_main.m

% Description: Main script for the Hoskins-West modified Eady-type model in the
% quasi-geostrophic framework. Initializes parameters, computes mean zonal wind
% (Ubar) and PV gradient (BPVy), constructs matrices for PV inversion and
% advection, solves the eigenvalue problem for wave stability, and saves results.

% This version uses the Hoskins–West Eady-type mean flow from
% Sample_basic_flows.pdf (p. 9):
%
%   U(y,z) = (g / (f0 * Theta0)) * (H * ΔT / Ly) *
%             [ (z/H) - (μ/2)*(z/H)^2
%               + sinh(2πLy/Lx * z/H)/sinh(2πLy/Lx) * cos(2π(y - y0)/Lx) ]
%
% where μ = 0 → Eady model; μ = 0.5 or 1 → Hoskins–West Eady-type.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx model
tic

addpath('functions'); % add functions folder
addpath('config');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load parameters from hwme_config.m
params = hwme_config();

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
y_s = params.y_s; % reference latitude (m)
N = params.N; % Brunt-Vaisala frequency (s^-1)
Lr = params.Lr; % Rossby radius of deformation (m)
Ly = params.Ly; % meridional domain length (50° latitude range, m)
prefac = params.prefac; % prefactor for Ubar calculation
hwme_data = fullfile(params.hwme_data_dir, params.hwme_data_filename);

%% Initialize Ubar and BPVy
Ubar = hwme_ubar(jj, kk, gg, f0, Theta0, dTbar, HH, Ly, ZZ, YY, mu, y_s, Lr);
BPVy = hwme_bpvy(jj, kk, beta, Ubar, dy, f0, NN2, dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% MATRIX CONSTRUCTION & EIGENVALUES %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct matrices B, C, D
[B, C, D] = hwme_matrices(ll, @stream2pv, @stream2xPVadv, @stream2yPVadv);

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

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(params.hwme_data_dir, 'dir')
    mkdir(params.hwme_data_dir);
end
save(hwme_data, 'Ubar', 'BPVy', 'eigVec', 'eigVal', 'eigVec2', 'eigVal2', ...
    'eigVec3', 'eigVal3', 'jj', 'kk', 'll', 'NN2', 'f0', 'dy', 'dz', 'm0', ...
    'Lx', 'beta', 'cplx', 'gg', 'Theta0', 'HH', 'B', 'C', 'D', 'model');

%% Plot Ubar
contourf(YY/1e6, ZZ/1e3, Ubar');
xlabel('y (×10^6 m)'); ylabel('z (km)');
title('Hoskins–West Eady-type Ū(y,z)');
colorbar
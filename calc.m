%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: calc.m
%
% DESCRIPTION: Main post-processing driver for the Hoskins-West modified Eady 
% (HWME) model. Loads eigenmode data, selects a mode, normalizes the 
% streamfunction, computes all diagnostic fields (geopotential, temperature, 
% winds, PV, vertical velocity), generates Hovmöller diagrams, and saves 
% results for plotting. No figures are generated.
%
% INPUT:
% - Eigenmode data: eigVec2, eigVal2, B, Ubar, BPVy (from main.m output)
% - Model parameters: from cfg()
%
% OUTPUT:
% - Saves structured results to: data/calc_wave-m0_eMode-n.mat
%
% MATH/FUNCTIONS:
% - phase_speed = -ω / kx,   growth_rate = σ,   e-folding = 1/σ (days)
% - Normalization: max|ψ| = 10 m²/s² (via gpt_h scaling)
% - All fields reconstructed via XV2field(·, ii, dx) * scaling
%
% VARIABLES:
% - XV: Streamfunction eigenvector (m²/s), length ll
% - omega: Imaginary part of eigenvalue → phase frequency (s⁻¹)
% - growth_rate: Real part of eigenvalue (s⁻¹)
% - phase_speed: -ω / (2π m0 / Lx) (m/s)
% - eFolding: e-folding time in days
% - QV: Potential vorticity anomaly = B * XV
% - XVy, XVx, XVz: Meridional, zonal, vertical derivatives of XV
% - gpt_h: Geopotential height anomaly (m)
% - temp: Temperature anomaly (K)
% - ug, vg: Geostrophic zonal and meridional wind (m/s)
% - pvfield: 3D PV anomaly field
% - w_vec: Interior vertical velocity from omega equation (m/s)
% - wfield: Full 3D vertical velocity field
% - *_hovmoler*: Time-longitude evolution at fixed latitude/level
% - Global params: jj, kk, ii, dx, dy, dz, f0, gg, Theta0, HH, m0, Lx, etc.
% - Functions used: XV2XVy, XV2XVx, XV2XVz, XV2field, Gmatrix, F1/F2/F3matrix, 
%   w2vec, w2field, XV2streamxtime, XV2ugxtime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll xx yy zz HH ii BPVy NN2 f0 gg dx dy m0 dz Lx Theta0 ...
    Ubar beta cplx hlat hlevel time LW n_mode fig_path

tic
addpath(['functions', filesep])

%% Load parameters
params = cfg();
jj = params.jj; kk = params.kk; ll = params.ll; Lx = params.Lx; LW = params.LW;
f0 = params.f0; beta = params.beta; gg = params.gg; Theta0 = params.Theta0;
HH = params.HH; dy = params.dy; dz = params.dz; NN2 = params.NN2;
cplx = params.cplx; m0 = params.m0; fig_path = params.plot_dir;
ii = params.ii; dx = params.dx; xx = params.xx; yy = params.yy; zz = params.zz;
hlat = params.hlat; time = params.time; n_mode = params.n_mode;

%% Load eigenmode data
data = fullfile(params.data_dir, params.data_filename);
load(data, 'eigVec2', 'eigVal2', 'B', 'Ubar', 'BPVy'); % Load necessary variables

fprintf('Loading computed data from %s for QG-CAI model...\nCalculating parameters (wfield goes brrrr)...\n', data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select mode
XV = eigVec2(:, n_mode);
omega = imag(eigVal2(n_mode));
phase_speed = -omega / (2 * pi * m0 / Lx);
growth_rate = real(eigVal2(n_mode));
eFolding = (1 / growth_rate) / 86400; % days

%% Normalise XV
valuemax = max(XV2field(XV, ii, dx) * f0 / gg, [], 'all');
XV = (10 / valuemax) * XV;

%% Compute fields
QV  = B * XV;
XVy = XV2XVy(XV);
XVx = XV2XVx(XV);
XVz = XV2XVz(XV);

gpt_h = XV2field(XV, ii, dx) * f0 / gg;
temp  = (f0 * HH / 287) * XV2field(XVz, ii, dx);
ug    = -XV2field(XVy, ii, dx);
vg    =  XV2field(XVx, ii, dx);
pvfield = XV2field(QV, ii, dx);

%% Create matrixes for vertical velocity
G  = Gmatrix();
F3 = F3matrix(XV);
F2 = F2matrix(XV,Ubar);
F1 = F1matrix(XV,Ubar);

w_vec   = w2vec(G,F1,F2,F3);
wfield  = w2field(w_vec, ii, dx);

%% Hovmoller diagrams calculations
gpt_h_hovmoler1  = XV2streamxtime(XV, ii, dx, omega, hlat, 1) * f0 / gg;
gpt_h_hovmoler51 = XV2streamxtime(XV, ii, dx, omega, hlat, 51) * f0 / gg;
ug_hovmoler1     = XV2ugxtime(XVy, ii, dx, omega, hlat, 1);
ug_hovmoler51    = XV2ugxtime(XVy, ii, dx, omega, hlat, 51);

XV = zeros(ll, 1);
XV(:) = eigVec2(:, n_mode);

%% Create output directory
if ~exist(params.plot_dir, 'dir')
    mkdir(params.plot_dir);
end

calcFile = fullfile(params.data_dir, params.calc_filename);

% ---- 1. Save -------------------------------------------------
save(calcFile, ...
    'gpt_h', 'temp', 'ug', 'vg', 'pvfield', 'wfield', ...
    'gpt_h_hovmoler1', 'gpt_h_hovmoler51', ...
    'ug_hovmoler1',  'ug_hovmoler51', ...
    'xx','yy','zz','time','Ubar','BPVy','XV', ...
    'm0','n_mode','growth_rate','omega','phase_speed','eFolding', ...
    'hlat','params');  

disp_str = ['All fields calculated and saved (individually) to ', calcFile];
disp(disp_str);
toc

save(fullfile(params.data_dir, params.calc_filename));

disp_str = ['All fields calculated and saved to ', params.data_dir, params.calc_filename];

disp(disp_str);
toc
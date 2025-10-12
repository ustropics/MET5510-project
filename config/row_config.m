%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: row_config.m

% Description: Defines constants and derived parameters for the Rossby wave 
% model in the quasi-geostrophic framework. Returns a structure containing 
% all parameters used for grid setup, physical constants, and model-specific 
% calculations.

% Variables defined in params structure:
% - cplx: Imaginary unit for complex number operations (sqrt(-1))
% - m0: Wave number (integer, typically 2 for planetary waves)
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - Lx: Zonal domain length at 45° latitude (m)
% - Ly: Meridional domain length (50° latitude range, m)
% - f0: Coriolis parameter at 45° latitude (s^-1)
% - beta: Beta parameter for meridional PV gradient (m^-1 s^-1)
% - gg: Gravitational acceleration (m s^-2)
% - Theta0: Reference potential temperature (K)
% - delta_Theta0: Potential temperature difference across domain (K)
% - HH: Scale height (m)
% - dy: Meridional grid spacing (m)
% - dz: Vertical grid spacing (m)
% - ll: Total number of linear indices for interior points
% - NN2: Brunt-Vaisala frequency squared (s^-2)
% - Ubar_const: Constant zonal wind speed (m s^-1)

% Directories and filenames:
% - row_data_dir: Directory for saving data
% - row_data_filename: Filename for saving data, includes wave number
% - row_plot_dir: Directory for saving plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = row_config()
params = struct();

%% Initial variables for complex and wave number
params.cplx = sqrt(-1); % imaginary unit for complex number operations
params.m0 = 2; % wave number (0-3 is longwave/planetary, 3+ is shortwave)

%% Grid point parameterization
params.jj = 50; % number of latitude grid points
params.kk = 50; % number of height grid points
params.Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)
params.Ly = 6.37 * 1.0e6 * 50 * pi/180; % meridional domain length (50° latitude range, m)

%% Coriolis, beta, and gravity calculations/constants
params.f0 = 2 * (2*pi/86400) * sin(pi/4); % coriolis parameter at 45° latitude (s^-1)
params.beta = 2 * (2 * pi/86400) * cos(pi/4)/(6.37*1.e6); % beta parameter (meridional PV gradient, m^-1 s^-1)
params.gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures, height, and additional domain gridding
params.Theta0 = 300; % reference potential temperature (K)
params.delta_Theta0 = 30; % potential temperature difference across domain (K)
params.HH = 10000; % scale height (m)

%% Derived grid parameters
params.dy = params.Ly / params.jj; % meridional grid spacing (m)
params.dz = params.HH / params.kk; % vertical grid spacing (m)
params.ll = (params.jj-1) * (params.kk+1); % total number of linear indices

%% Brunt-Vaisala frequency
params.NN2 = params.gg * params.delta_Theta0 / (params.HH * params.Theta0); % s^-2

%% Constant zonal wind
params.Ubar_const = 10; % constant zonal wind speed (m/s)

%% Directory variables
params.row_plot_dir = 'output/plots/'; % directory for saving plots
params.row_data_dir = 'output/data/'; % directory for saving data
params.row_data_filename = ['row_wave_', num2str(params.m0), '.mat']; % data filename

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: eady_config.m

% Description: Defines constants and derived parameters for the Eady wave 
% model in the quasi-geostrophic framework. Returns a structure containing 
% all parameters used for grid setup, physical constants, and plotting.

% Variables defined in params structure:
% - cplx: Imaginary unit for complex number operations
% - m0: Wave number
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - ll: Total number of linear indices
% - Lx: Zonal domain length (m)
% - Ly: Meridional domain length (m)
% - f0: Coriolis parameter (s^-1)
% - beta: Beta parameter (m^-1 s^-1, set to 0 for Eady model)
% - gg: Gravitational acceleration (m s^-2)
% - Theta0: Reference potential temperature (K)
% - delta_Theta0: Potential temperature difference (K)
% - HH: Scale height (m)
% - dy: Meridional grid spacing (m)
% - dz: Vertical grid spacing (m)
% - NN2: Brunt-Vaisala frequency squared (s^-2)
% - ZZ: Vertical grid coordinates (m)
% - YY: Meridional grid coordinates (m)
% - dTbar: Potential temperature gradient (K)
% - Ubar0: Base zonal wind speed (m/s)
% - ii: Number of longitude grid points
% - dx: Longitudinal grid spacing (m)
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - zz: Height coordinates (km)
% - hlat: Latitude index for Hovmoller diagram
% - hlevel: Vertical level index for Hovmoller diagram
% - time: Time coordinates for Hovmoller (days)
% - n_mode: Mode number for plotting
% - eady_data_dir: Directory for saving data
% - eady_data_filename: Filename for saving data
% - eady_plot_dir: Directory for saving plots
% - Lambda: Vertical shear (s^-1)
% - qy_surf: Surface PV gradient (s^-1 m^-1)
% - qy_trop: Tropopause PV gradient (s^-1 m^-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = eady_config()
params = struct();

%% Initial variables for complex and wave number
params.model = 'eady';
params.cplx = sqrt(-1); % imaginary unit
params.m0 = 2; % wave number

%% Grid point parameterization
params.jj = 50; % number of latitude grid points
params.kk = 50; % number of height grid points
params.Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)
params.Ly = 6.37 * 1.0e6 * 50 * pi/180; % meridional domain length (50° latitude range, m)


%% Coriolis, beta, and gravity constants
params.f0 = 2 * (2 * pi / 86400) * sin(pi / 4); % Coriolis parameter at 45° latitude (s^-1)
params.beta = 0; % beta parameter (m^-1 s^-1, set to 0 for Eady model)
params.gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures and height
params.Theta0 = 300; % reference potential temperature (K)
params.delta_Theta0 = 60; % potential temperature difference (K)
params.HH = 10000; % scale height (m)

%% Derived grid parameters
params.dy = params.Ly / params.jj; % meridional grid spacing (m)
params.dz = params.HH / params.kk; % vertical grid spacing (m)
params.ll = (params.jj - 1) * (params.kk + 1); % total number of linear indices

%% Brunt-Vaisala frequency
params.NN2 = params.gg * params.delta_Theta0 / (params.HH * params.Theta0); % s^-2

%% Grid coordinates
params.ZZ = 0.0:params.dz:params.HH; % vertical grid (m)
params.YY = 6.37 * 1.0e6 * pi / 4 - params.Ly / 2 : params.dy : 6.37 * 1.0e6 * pi / 4 + params.Ly / 2; % meridional grid (m)

%% Constants for Ubar initialization
params.dTbar = 60; % potential temperature gradient (K)
params.Ubar0 = 0; % base zonal wind speed (m/s)

%% Background flow parameters
params.Lambda = (params.gg / params.f0 / params.Theta0) * (params.dTbar / params.Ly); % vertical shear (s^-1)
params.qy_surf = - (params.f0^2 / (params.NN2 * params.HH)) * params.Lambda; % surface PV gradient (s^-1 m^-1)
params.qy_trop = (params.f0^2 / (params.NN2 * params.HH)) * params.Lambda; % tropopause PV gradient (s^-1 m^-1)

%% Plotting parameters
params.ii = 360; % longitude grid
params.dx = params.Lx / params.ii; % longitudinal grid spacing (m)
params.xx = 0.0 : 360 / params.ii : 360; % longitude coordinates (degrees)
params.yy = linspace(45 - 25, 45 + 25, params.jj + 1); % latitude coordinates (degrees)
params.zz = linspace(0.0, 10, params.kk + 1); % height coordinates (km)
params.hlat = floor(params.jj / 4 + 1); % lat for Hovmoller diagram
params.hlevel = 1; % vertical level for Hovmoller diagram
params.time = 0 : 1 : 50; % time coordinates (days)
params.n_mode = 7; % mode number

%% Directory variables
params.eady_plot_dir = 'output/figures/'; % directory for saving plots
params.eady_data_dir = 'output/data/'; % directory for saving data
params.eady_data_filename = ['eady_wave_', num2str(params.m0), '.mat']; % data filename

end
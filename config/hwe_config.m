%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: hwe_config.m

% Description: Defines constants and derived parameters for the modified 
% Hoskins-West Eady-type Model in the quasi-geostrophic framework. Returns a 
% structure containing all parameters used for grid setup, physical constants, 
% and model-specific calculations.

% Variables defined in params structure:
% - cplx: Imaginary unit for complex number operations (sqrt(-1))
% - m0: Wave number (integer, typically 7)
% - jj: Number of latitude grid points
% - kk: Number of height grid points
% - Lx: Zonal domain length at 45° latitude (m)
% - Ly: Meridional domain length (50° latitude range, m)
% - f0: Coriolis parameter at 45° latitude (s^-1)
% - beta: Beta parameter for meridional PV gradient (m^-1 s^-1)
% - gg: Gravitational acceleration (m s^-2)
% - Theta0: Reference potential temperature (K)
% - DeltaTheta: Vertical/horizontal potential temperature difference (K)
% - HH: Scale height (m)
% - dy: Meridional grid spacing (m)
% - dz: Vertical grid spacing (m)
% - ll: Total number of linear indices for interior points
% - NN2: Brunt-Vaisala frequency squared (s^-2)
% - DeltaT_hor: Horizontal temperature difference (K)
% - mu: Mu parameter for quadratic term in zonal wind (unused in Eady)
% - U0: Mean zonal wind offset (m s^-1)
% - L_d: Rossby deformation radius
% - m_y: Meridional wavenumber (unused in Eady)
% - gamma: Gamma parameter for sinh term (unused in Eady)
% - y0: Meridional domain center (m)
% - Lambda: Shear rate for zonal wind calculation (derived from DeltaT_hor)

% Directories and filenames:
% - hwe_plot_dir: Directory for saving plots
% - hwe_data_dir: Directory for saving data
% - hwe_data_filename: Filename for saving data, includes wave number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = hwe_config()
params = struct();

%% Initial variables for complex and wave number
params.cplx = sqrt(-1); % imaginary unit for complex number operations
params.m0 = 1; % wave number
params.n_mode = 7; 

%% Grid point parameterization
params.jj = 50; % number of latitude grid points
params.kk = 50; % number of height grid points
params.Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)
params.Ly = 6.37 * 1.0e6 * 50 * pi/180; % meridional domain length (50° latitude range, m)

%% Coriolis, beta, and gravity calculations/constants
params.f0 = 2 * (2*pi/86400) * sin(pi/4); % coriolis parameter at 45° latitude (s^-1)
params.beta = 0; % beta parameter set to 0 for f-plane Eady-type (m^-1 s^-1)
params.gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures, height, and additional domain gridding
params.Theta0 = 300; % reference potential temperature (K)
params.DeltaTheta = 40; % vertical/horizontal potential temperature difference (K), adjustable to 30 or 40
params.HH = 10000; % scale height (m)

%% Derived grid parameters
params.dy = params.Ly / params.jj; % meridional grid spacing (m)
params.dz = params.HH / params.kk; % vertical grid spacing (m)
params.ll = (params.jj-1) * (params.kk+1); % total number of linear indices

%% Brunt-Vaisala frequency
params.NN2 = params.gg * params.DeltaTheta / (params.HH * params.Theta0); % s^-2, matching N^2_* from slide

%% Parameters for Modified Hoskins-West Eady-type Model
params.DeltaT_hor = params.DeltaTheta; % horizontal temperature difference (K), set to DeltaTheta for consistency
params.mu = 0; % mu parameter (set to 0 as per slide for Eady, unused in new ubarCalc)
params.U0 = 5; % mean zonal wind offset (m/s) to match previous context and slide flexibility
params.L_d = sqrt(params.NN2 * params.HH / params.f0); % Rossby deformation radius, matching L_R from slide
params.m_y = 2 * pi / params.Ly; % meridional wavenumber (kept for compatibility, though unused)
params.gamma = params.m_y * params.L_d; % gamma for sinh term (kept but unused)
params.y0 = params.Ly / 2; % meridional domain center (m), assuming y_s = Ly/2 for simplicity

%% Directory variables
params.hwe_plot_dir = 'output/plots/'; % directory for saving plots
params.hwe_data_dir = 'output/data/'; % directory for saving data
params.hwe_data_filename = ['hwe_wave_', num2str(params.m0), '.mat']; % data filename

end
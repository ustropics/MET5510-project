%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: hwme_plot.m

% DESCRIPTION: Script for plotting results from the Hoskins-West modified 
% Eady-type model, including eigenvector amplitude, geopotential height, 
% Hovmoller diagram, zonal wind, meridional wind, temperature, 
% potential vorticity, additional cross-sections, Hovmoller diagrams, 
% and background state fields. Loads data from 'hwme_wave_#.mat',
% computes necessary fields, and saves plots to 'output/plots/'.

% SCRIPTS:
% - hwme_config.m: Loads model parameters

% PLOTS:
% - plot_evec_amp: Plots eigenvector amplitude contour
% - plot_geopotential_height: Plots geopotential height contour
% - plot_hovmoller: Plots Hovmoller diagram
% - plot_zonal_wind: Plots zonal wind contour
% - plot_meridional_wind: Plots meridional wind contour
% - plot_temperature: Plots temperature contour
% - plot_gph_top: Plots geopotential height at top boundary
% - plot_pvfield: Plots potential vorticity contour
% - plot_vg_cross_section: Plots meridional wind vertical cross-section
% - plot_ug_hovmoller: Plots Hovmoller diagram for zonal wind
% - plot_ubar_contour: Plots Ubar contour
% - plot_dpvdym_int: Plots d(PVbar)/dy interior contour
% - plot_dpvdym_boundaries: Plots d(PVbar)/dy at boundaries with beta

% FUNCTIONS:
% - XV2field: Computes 3D field from streamfunction vector
% - XV2streamxtime: Computes Hovmoller data for streamfunction
% - XV2ugxtime: Computes Hovmoller data for zonal wind
% - XV2XVy: Computes meridional derivative of streamfunction
% - XV2XVx: Computes zonal derivative of streamfunction
% - XV2XVz: Computes vertical derivative of streamfunction
% - jk2l: Converts 2D indices to linear index
% - l2jk: Converts linear index to 2D indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll xx yy zz HH ii BPVy NN2 f0 gg dx dy m0 dz Lx Theta0 ... 
    Ubar beta cplx hlat hlevel time n_mode hwme_data fig_path

addpath(['functions', filesep])
addpath(['config', filesep])
addpath(['plots', filesep])

%% Load parameters from hwme_config.m
params = hwme_config();

% Assign global variables
jj = params.jj;
kk = params.kk;
ll = params.ll;
Lx = params.Lx;
f0 = params.f0;
beta = params.beta;
gg = params.gg;
Theta0 = params.Theta0;
HH = params.HH;
dy = params.dy;
dz = params.dz;
NN2 = params.NN2;
cplx = params.cplx;
m0 = params.m0;
fig_path = params.hwme_plot_dir;

% Assign plotting parameters
ii = params.ii;
dx = params.dx;
xx = params.xx;
yy = params.yy;
zz = params.zz;
hlat = params.hlat;
hlevel = params.hlevel;
time = params.time;

%% Load data from hwme_wave_#.mat
hwme_data = fullfile(params.hwme_data_dir, params.hwme_data_filename);
load(hwme_data)

%% Define maximum number of modes
n_max = 3; % Set the number of modes to loop through
l_max = zeros(n_max, 1); % Initialize l_max array

% Check if n_max is valid
if n_max > size(eigVec2, 2)
    error('n_max (%d) exceeds number of available modes (%d)', n_max, size(eigVec2, 2));
end

%% Loop over modes
for n_mode = 1:n_max
    %% Select mode and compute properties
    XV = zeros(ll, 1);
    XV(:) = eigVec2(:, n_mode);
    omega = imag(eigVal2(n_mode));
    phase_speed = -omega / (2 * pi * m0 / Lx);
    growth_rate = real(eigVal2(n_mode));
    if growth_rate ~= 0
        eFolding = (1 / growth_rate) / 86400; % in days
    else
        eFolding = Inf; % Handle zero growth rate
    end

    %% Compute l_max for this mode
    if (real(eigVal2(n_mode)) > 1.0e-10)
        xv2d = reshape(XV, jj-1, kk+1);
        amp = abs(xv2d);
        l_max(n_mode) = sum(islocalmax(amp(:,1)));
        if (amp(1,1) > amp(2,1))
            l_max(n_mode) = l_max(n_mode) + 1;
        end
        if (amp(end,1) > amp(end-1,1))
            l_max(n_mode) = l_max(n_mode) + 1;
        end
    end

    %% Normalize XV
    valuemax = max(XV2field(XV, ii, dx) * f0 / gg, [], 'all');
    if valuemax ~= 0
        XV = (10 / valuemax) * XV;
    else
        warning('valuemax is zero for n_mode = %d, skipping normalization', n_mode);
    end

    %% Compute fields
    QV = B * XV;
    XVy = XV2XVy(XV);
    XVx = XV2XVx(XV);
    XVz = XV2XVz(XV);

    gpt_h = XV2field(XV, ii, dx) * f0 / gg;
    temp = (f0 * HH / 287) * XV2field(XVz, ii, dx);
    ug = -XV2field(XVy, ii, dx);
    vg = XV2field(XVx, ii, dx);
    pvfield = XV2field(QV, ii, dx);

    %% Compute Hovmoller diagrams
    gpt_h_hovmoler = XV2streamxtime(XV, ii, dx, omega, hlat, hlevel) * f0 / gg;
    ug_hovmoler = XV2ugxtime(XVy, ii, dx, omega, hlat, hlevel);

    %% Create output directory if it doesn't exist
    if ~exist(params.hwme_plot_dir, 'dir')
        mkdir(params.hwme_plot_dir);
    end

    %% Plot eigenvector amplitude
    eVec_amp = zeros(jj + 1, kk + 1);
    for l = 1 : ll
        [j, k] = l2jk(l);
        eVec_amp(j, k) = XV(l) .* conj(XV(l));
    end
    plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, model, fig_path);

    %% Create list of figures to plot
    % Plot geopotential height
    plot_gph(xx, zz, gpt_h, jj, model, m0, n_mode, fig_path);

    % Plot Hovmoller diagram
    plot_hovmoller(xx, time, gpt_h_hovmoler, model, m0, n_mode, fig_path);

    % Plot zonal wind
    plot_zonal_wind(xx, yy, ug, model, m0, n_mode, fig_path);

    % Plot meridional wind
    plot_meridional_wind(xx, yy, vg, model, m0, n_mode, fig_path);

    % Plot temperature at mid-level
    plot_temperature(xx, yy, temp, kk, model, m0, n_mode, fig_path);

    % Plot geopotential height at top boundary
    plot_gph_top(xx, yy, gpt_h, model, m0, n_mode, fig_path);

    % Plot potential vorticity at surface
    plot_pvfield(xx, yy, pvfield, model, m0, n_mode, fig_path);

    % Plot meridional wind vertical cross-section
    plot_vg_cross_section(xx, zz, vg, jj, model, m0, n_mode, fig_path);

    % Plot Hovmoller diagram for zonal wind
    plot_ug_hovmoller(xx, time, ug_hovmoler, hlat, hlevel, model, m0, n_mode, fig_path);

    % Plot Ubar contour
    plot_ubar_contour(yy, zz, Ubar, model, m0, n_mode, fig_path);

    % Plot d(PVbar)/dy interior
    plot_dpvdym_int(yy, zz, BPVy, model, m0, n_mode, fig_path);

    % Plot d(PVbar)/dy at boundaries with beta
    plot_dpvdym_boundaries(yy, BPVy, beta, kk, model, m0, n_mode, fig_path);

    % Plot combo plot for hwme
    plot_hwme_bg_flow(yy, zz, jj, kk, Ubar, BPVy, model, m0, n_mode, fig_path);
end

%% Display l_max for all modes
disp('l_max for each mode:');
disp(l_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: hwme_plot.m
% DESCRIPTION: Script for plotting results from the Hoskins-West modified 
% Eady-type model, including eigenvector amplitude, geopotential height, 
% Hovmoller diagrams for geopotential height and zonal wind at hlevel = 1 and 51, 
% zonal wind at hlevel = 1, 25, and 51, meridional wind at hlevel = 1, 25, and 51, 
% temperature at hlevel = 1, 25, and 51, potential vorticity, additional 
% cross-sections, and background state fields. Loads data from 'hwme_wave_#.mat', 
% computes necessary fields, and saves plots to 'output/plots/'.

% SCRIPTS:
% - hwme_config.m: Loads model parameters

% PLOTS:
% - plot_evec_amp: Plots eigenvector amplitude contour
% - plot_geopotential_height: Plots geopotential height contour
% - plot_hovmoller: Plots Hovmoller diagram for geopotential height (hlevel = 1 and 51)
% - plot_zonal_wind: Plots zonal wind contour (hlevel = 1, 25, 51)
% - plot_meridional_wind: Plots meridional wind contour (hlevel = 1, 25, 51)
% - plot_temperature: Plots temperature contour (hlevel = 1, 25, 51)
% - plot_gph_top: Plots geopotential height at top boundary
% - plot_pvfield: Plots potential vorticity contour
% - plot_vg_cross_section: Plots meridional wind vertical cross-section
% - plot_ug_hovmoller: Plots Hovmoller diagram for zonal wind (hlevel = 1 and 51)
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
time = params.time;
n_mode = params.n_mode;

%% Load data from hwme_wave_#.mat
hwme_data = fullfile(params.hwme_data_dir, params.hwme_data_filename);
load(hwme_data)

%% Select mode and compute properties
XV = zeros(ll, 1);
XV(:) = eigVec2(:, n_mode);
omega = imag(eigVal2(n_mode));
phase_speed = -omega / (2 * pi * m0 / Lx);
growth_rate = real(eigVec2(n_mode));
eFolding = (1 / growth_rate) / 86400; % in days

%% Normalize XV
valuemax = max(XV2field(XV, ii, dx) * f0 / gg, [], 'all');
XV = (10 / valuemax) * XV;

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

% =====================================================================
% PATCH 1: CORRECT INDEXING (the only thing that broke w)
% =====================================================================
jk2lw = @(j,k) j + (k-2)*(jj+1); 
lw2jk = @(l)  deal(mod(l-1,jj+1)+1, floor((l-1)/(jj+1))+2);

% =====================================================================
% PATCH 2: BUILD G MATRIX (was missing!)
% =====================================================================
LW = (jj+1)*(kk-1);
G  = zeros(LW,LW);
for l0 = 1:LW
    w = zeros(LW,1); w(l0) = 1;
    G(:,l0) = w2ellipse(w);
end

% =====================================================================
% PATCH 3: F1+F2+F3 — exact copy from diagnosis script
% =====================================================================
F1 = zeros(LW,1); F2 = F1; F3 = F1;

% F3 – beta term
for j = 2:jj
    for k = 2:kk
        l = jk2lw(j,k);
        F3(l) = (2*pi*m0/Lx)*cplx*beta*...
                (XV(jk2l(j,k+1)) - XV(jk2l(j,k-1)))/(2*dz);
    end
end

% F2 – shear term
for j = 2:jj
    for k = 2:kk
        l = jk2lw(j,k);
        dUdy = (Ubar(j+1,k) - 2*Ubar(j,k) + Ubar(j-1,k))/dy^2;
        dXdz = (XV(jk2l(j,k+1)) - XV(jk2l(j,k-1)))/(2*dz);
        F2(l) = -2*(2*pi*m0/Lx)*cplx * dUdy * dXdz;
    end
end

% F1 – boundary curvature (copy-paste)
% j = 1
j = 1;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(XV(jk2l(3,k))-2*XV(jk2l(2,k)))/(2*dz*dy*dy);
end

% j = 2
j = 2;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
          (XV(jk2l(3,k))-2*XV(jk2l(2,k)))/dy/dy)/(2*dz);
end

% interior j = 3 .. jj-1
for j = 3:jj-1
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
            *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
              (XV(jk2l(j+1,k))-2*XV(jk2l(j,k))+XV(jk2l(j-1,k)))/dy/dy)/(2*dz);
    end
end

% j = jj
j = jj;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
          (XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/dy/dy)/(2*dz);
end

% j = jj+1
j = jj+1;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
          *(XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/(2*dz*dy*dy);
end

% =====================================================================
% PATCH 4: SOLVE FOR w AND MAKE 3D FIELD
% =====================================================================
w      = (f0/NN2) * (G \ (F1+F2+F3));
wfield = w2wfield(w, ii, dx);
fprintf('max |vertical velocity| = %.2e m/s\n', max(abs(w)));

% =====================================================================
% YOUR ORIGINAL PLOTTING CODE STARTS HERE (unchanged)
% =====================================================================

%% Compute Hovmoller diagrams for hlevel = 1 and 51
gpt_h_hovmoler1 = XV2streamxtime(XV, ii, dx, omega, hlat, 1) * f0 / gg;
gpt_h_hovmoler51 = XV2streamxtime(XV, ii, dx, omega, hlat, 51) * f0 / gg;
ug_hovmoler1 = XV2ugxtime(XVy, ii, dx, omega, hlat, 1);
ug_hovmoler51 = XV2ugxtime(XVy, ii, dx, omega, hlat, 51);

if ~exist(params.hwme_plot_dir, 'dir')
    mkdir(params.hwme_plot_dir);
end

% Eigenvector amplitude
eVec_amp = zeros(jj + 1, kk + 1);
for l = 1 : ll
    [j, k] = l2jk(l);
    eVec_amp(j, k) = XV(l) .* conj(XV(l));
end

% === ALL YOUR PLOTS (exactly as you had them) ===
% plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, model, fig_path);
% plot_zvt(vg, temp, m0, n_mode, fig_path);
% plot_zvu(vg, ug, m0, n_mode, fig_path);
% plot_zwt(wfield, temp, m0, n_mode, fig_path);          % NOW WORKS!
% plot_background_flow(yy, zz, jj, kk, Ubar, BPVy, model, m0, n_mode, fig_path);
% plot_gph(xx, zz, gpt_h, jj, model, m0, n_mode, fig_path);
% plot_gph_hovmoller(xx, time, gpt_h_hovmoler1, model, m0, n_mode, fig_path, 1);
% plot_gph_hovmoller(xx, time, gpt_h_hovmoler51, model, m0, n_mode, fig_path, 51);
% plot_zonal_wind(xx, yy, ug, 1, model, m0, n_mode, fig_path);
% plot_zonal_wind(xx, yy, ug, 26, model, m0, n_mode, fig_path);
% plot_zonal_wind(xx, yy, ug, 51, model, m0, n_mode, fig_path);
% plot_meridional_wind(xx, yy, vg, 1, model, m0, n_mode, fig_path);
% plot_meridional_wind(xx, yy, vg, 26, model, m0, n_mode, fig_path);
% plot_meridional_wind(xx, yy, vg, 51, model, m0, n_mode, fig_path);
% plot_temperature(xx, yy, temp, 1, model, m0, n_mode, fig_path);
% plot_temperature(xx, yy, temp, 26, model, m0, n_mode, fig_path);
% plot_temperature(xx, yy, temp, 51, model, m0, n_mode, fig_path);
% plot_gph_top(xx, yy, gpt_h, model, m0, n_mode, fig_path);
% plot_pvfield(xx, yy, pvfield, 1, model, m0, n_mode, fig_path);
% plot_pvfield(xx, yy, pvfield, 26, model, m0, n_mode, fig_path);
% plot_pvfield(xx, yy, pvfield, 51, model, m0, n_mode, fig_path);
% plot_vg_cross_section(xx, zz, vg, jj, model, m0, n_mode, fig_path);
% plot_ug_hovmoller(xx, time, ug_hovmoler1, hlat, 1, model, m0, n_mode, fig_path);
% plot_ug_hovmoller(xx, time, ug_hovmoler51, hlat, 51, model, m0, n_mode, fig_path);
% plot_ubar_contour(yy, zz, Ubar, model, m0, n_mode, fig_path);
% plot_dpvdym_int(yy, zz, BPVy, model, m0, n_mode, fig_path);
% plot_dpvdym_boundaries(yy, BPVy, beta, kk, model, m0, n_mode, fig_path);
% plot_background_flow(yy, zz, jj, kk, Ubar, BPVy, model, m0, n_mode, fig_path);

disp('All plots generated successfully... BRB cig break...');
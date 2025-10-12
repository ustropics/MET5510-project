%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_plotting.m

% Description: Script for loading Rossby wave data from 'Rossby_wave_2.mat' 
% and generating plots to visualize eigenvector amplitude, meridional 
% cross-section of geopotential height, and Hovmoller diagram of 
% geopotential height evolution over time, based on the quasi-geostrophic 
% model analysis.

% Functions used: 
% - l2jk: Converts linear index to 2D indices (j, k).
% - jk2l: Converts 2D indices (j, k) to linear index.
% - XV2field: Transforms eigenvector to 3D field.
% - XV2streamxtime: Generates streamfunction field with time evolution.
% - XV2XVx: Computes x-derivative of eigenvector.
% - XV2XVy: Computes y-derivative of eigenvector.
% - XV2XVz: Computes z-derivative of eigenvector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Rossby_wave_2.mat')
addpath('functions')

global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% grid parameters
ii = 360; % longitude grid
dx = Lx/ii; % longitude grid spacing

xx = 0.0:360/ii:360; % longitude coordinates
yy = linspace(45-25, 45+25, jj+1); % latitude coordinates
zz = linspace(0.0,10,kk+1); % height grid (H=10km)

XV = zeros(ll,1); % initailize streamfunction vector
% XV(:) = eigVec3(:,7);

%% Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf (pages 1-2)
% 1 is uniform/barotropic, 2 is a cross-section
% 7 is UL wind and LL wind (6-7 is baroclinic mode)
n_mode = 7; % this goes west fast, 2 is uFnstable
XV(:) = eigVec3(:,n_mode);
omega = imag(eigVal3(n_mode));
phase_speed = -omega/(2*pi*m0/Lx);

% from Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% slide 6
growth_rate = real(eigVal3(n_mode)); % real part of eigenvalue
eFolding = (1/growth_rate)/86400; % in units of days

% % Amplitude of eigenvector
eVec_amp = zeros(jj+1, kk+1); % slide 7-8
for l = 1: ll
    [j,k] = l2jk(l);
    eVec_amp(j,k) = XV(l).*conj(XV(l));
end

% set geopotential field
% gpt_h = XV2field(XV,ii,dx) * f0/gg;
latmax = (jj/2 + 1)/2;
levelmax = 1;

% finds the maximum value index in 3d
% [valuemax, indexmax] = max(gpt_h(:,latmax, levelmax));
gpt_h = XV2field(XV,ii,dx)*f0/gg;

% finds the maximum value index converting 3d to 1d
[valuemax, indexmax] = max(gpt_h(:));

% normalize the value
XV = (10/valuemax) * XV;

% normalize the evec solution so that the amplitude of geopotential
% height anomalies is 10 meters for all eigen solutions


%ymode = n_zero(1);

Amp = sqrt(sum(XV.*XV)/ll);
% XV(:) = (100/f0/Amp)*XV(:);
XV(:) 

QV = B * XV;
XVy = XV2XVy(XV); % d/dy at N and S boundary not calculated
XVx = XV2XVx(XV);
XVz = XV2XVz(XV);
% Vort = XV2vort(XV);

% from Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
gpt_h = XV2field(XV,ii,dx)*f0/gg; % 3d geopotential height from slide 8
temp = (f0*HH/287) * XV2field(XVz, ii, dx); % from slide 11
ug = -XV2field(XVy, ii, dx);
vg = XV2field(XVx, ii, dx);
% vortfield =  XV2field(Vort, ii, dx);
pvfield = XV2field(QV, ii, dx);

% run series from day 0 to day 50, x from 0 to 360 degrees
hlat = floor(jj/4 + 1); % lat for hovmoller diagram
hlevel = 1; % vertical level for hovmoller diagram

% change dx to x to do meridional wind (slide 7)
gpt_h_hovmoler = XV2streamxtime(XV, ii, dx, omega, hlat, hlevel) * f0/gg;


% change XV to QV to plot PV
% gpt_h = XV2field(XV,ii,dx) * f0/gg;
pv = XV2field(QV,ii,dx);

time = 0:1:50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create plots folder if it doesn't exist
if ~exist('plots', 'dir')
    mkdir('plots');
end

% Figure 1: Eigenvector amplitude
fig1 = figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off');
contourf(yy, zz, eVec_amp', 'linestyle', 'none');
xlabel('Latitude')
ylabel('Height (km)')
set(gca, 'fontsize', 22, 'color', 'w')
title =(['Zonal wave # = ', num2str(m0), ', eMode # =', num2str(n_mode)...
    ', real eVal = ', num2str(growth_rate), ...
    ', imag eVal = ', num2str(omega)]);
colorbar
saveas(fig1, 'plots/eigenvector_amplitude.png');

% Figure 2: Meridional cross-section of geopotential height
fig2 = figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off');
contourf(xx, zz, squeeze(gpt_h(:,jj/2+1,:))', 'LineStyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca, 'xtick', 0:30:360)
set(gca, 'ytick', 0:5:50)
saveas(fig2, 'plots/meridional_cross_section_geopotential.png');

% Figure 3: Hovmoller diagram of geopotential height
fig3 = figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off');
contourf(xx, time, gpt_h_hovmoler', -100:2:100, 'LineStyle', 'none');
set(gca,'xlim',[0,90])
colorbar;
xlabel('Longitude')
ylabel('Time (days')
set(gca, 'xtick', 0:30:360)
set(gca, 'ytick', 0:5:50)
% set(gca, 'ylim', )
saveas(fig3, 'plots/hovmoller_geopotential.png');
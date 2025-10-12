%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: row_plot.m

% Description: Script for loading Rossby wave data from 'row_wave_#.mat' 
% and generating plots to visualize eigenvector amplitude, meridional 
% cross-section of geopotential height, Hovmoller diagram of geopotential 
% height evolution, and vertical motion diagnostics including geopotential 
% height, temperature, wind components, and vertical velocity fields.

% Functions used: 
% - l2jk: Converts linear index to 2D indices (j, k).
% - jk2l: Converts 2D indices (j, k) to linear index.
% - jk2lw: Converts 2D indices (j, k) to a weighted linear index.
% - XV2field: Transforms eigenvector to 3D field.
% - XV2streamxtime: Generates streamfunction field with time evolution.
% - XV2XVx: Computes x-derivative of eigenvector.
% - XV2XVy: Computes y-derivative of eigenvector.
% - XV2XVz: Computes z-derivative of eigenvector.
% - XVz2field: Converts z-derivative of eigenvector to 3D field.
% - XVy2field: Converts y-derivative of eigenvector to 3D field.
% - XVx2field: Converts x-derivative of eigenvector to 3D field.
% - w2ellipse: Transforms weight vector to an elliptical representation.
% - w2wfield: Converts weight vector to a 3D wind field.
% - compute_F123: Computes forcing terms F1, F2, F3 for vertical motion.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('functions'); % add functions folder
addpath('.'); % add current directory for row_config.m

% Load constants from row_config.m for directories and parameters
params = row_config();

global jj kk ll LW BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data from row_wave_#.mat
row_data_dir = params.row_data_dir;
row_plot_dir = params.row_plot_dir;
row_data_filename = params.row_data_filename;
row_data = fullfile(row_data_dir, row_data_filename);
load(row_data);

%% Grid parameters
ii = 360; % longitude grid points
dx = Lx/ii; % longitude grid spacing

xx = 0.0:360/ii:360; % longitude coordinates
yy = linspace(45-25, 45+25, jj+1); % latitude coordinates
zz = linspace(0.0,10,kk+1); % height grid (H=10km)

%% For vertical motion
LW = (jj+1)*(kk-1); % matrix size for vertical motion

%% Select mode
n_mode = 7; % example mode, adjust as needed

%% Initialize XV and compute basic fields
XV = zeros(ll,1); % initialize streamfunction vector
XV(:) = eigVec3(:,n_mode);
omega = imag(eigVal3(n_mode));
phase_speed = -omega/(2*pi*m0/Lx);
growth_rate = real(eigVal3(n_mode)); % real part of eigenvalue
eFolding = (1/growth_rate)/86400; % in units of days

%% Amplitude of eigenvector
eVec_amp = zeros(jj+1, kk+1);
for l = 1: ll
    [j,k] = l2jk(l);
    eVec_amp(j,k) = XV(l).*conj(XV(l));
end

%% Normalize XV for visualization
gpt_h = XV2field(XV,ii,dx) * f0/gg;
[valuemax, ~] = max(abs(gpt_h(:)));
XV = (10/valuemax) * XV;

%% Recompute fields after normalization
QV = B * XV;
XVy = XV2XVy(XV);
XVx = XV2XVx(XV);
XVz = XV2XVz(XV);

gpt_h = XV2field(XV,ii,dx)*f0/gg; % geopotential height
temp = (f0*HH/287) * XV2field(XVz, ii, dx); % temperature
ug = -XV2field(XVy, ii, dx); % zonal wind
vg = XV2field(XVx, ii, dx); % meridional wind
pvfield = XV2field(QV, ii, dx); % PV field

%% Hovmoller diagram setup
hlat = floor(jj/4 + 1); % lat for hovmoller diagram
hlevel = 1; % vertical level for hovmoller diagram
gpt_h_hovmoler = XV2streamxtime(XV, ii, dx, omega, hlat, hlevel) * f0/gg;
time = 0:1:50;

%% Vertical motion computation
% Initialize G matrix for vertical motion
G=zeros(LW,LW);
for l0 = 1:LW
    w=zeros(LW,1);
    w(l0)=1;
    
    % G matrix
    EW = w2ellipse(w);
    G(:,l0)=EW(:);
end

%% Compute F1, F2, and F3 using dedicated function
[F1, F2, F3] = row_F123(params, XV, Ubar);

%% Vertical velocity (w) using the elliptic equation: G*w = F1 + F2 + F3
w=(f0/NN2)*(G^-1)*(F1+F2+F3);
wfield=w2wfield(w,ii,dx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create plots folder if it doesn't exist
if ~exist(row_plot_dir, 'dir')
    mkdir(row_plot_dir);
end

%% Plot: Eigenvector amplitude
fig1 = figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off');
contourf(yy, zz, eVec_amp', 'linestyle', 'none');
xlabel('Latitude')
ylabel('Height (km)')
set(gca, 'fontsize', 22, 'color', 'w')
title_str =(['Zonal wave # = ', num2str(m0), ', eMode # =', num2str(n_mode)...
    ', real eVal = ', num2str(growth_rate), ...
    ', imag eVal = ', num2str(omega)]);
title(title_str);
colorbar
ROW_EVEC_AMP = fullfile(row_plot_dir, 'row_eigenvector_amplitude.png');
saveas(fig1, ROW_EVEC_AMP);

%% Plot: Meridional cross-section of geopotential height
fig2 = figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off');
contourf(xx, zz, squeeze(gpt_h(:,jj/2+1,:))', 'LineStyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca, 'xtick', 0:30:360)
set(gca, 'ytick', 0:5:50)
ROW_MERID_GPH = fullfile(row_plot_dir, 'row_meridional_cross_section_geopotential.png');
saveas(fig2, ROW_MERID_GPH);

%% Plot: Hovmoller diagram of geopotential height
fig3 = figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off');
contourf(xx, time, gpt_h_hovmoler', -100:2:100, 'LineStyle', 'none');
set(gca,'xlim',[0,90])
colorbar;
xlabel('Longitude')
ylabel('Time (days)')
set(gca, 'xtick', 0:30:360)
set(gca, 'ytick', 0:5:50)
ROW_HOVMOLLER = fullfile(row_plot_dir, 'row_hovmoller_geopotential.png');
saveas(fig3, ROW_HOVMOLLER);

%% Combined Plot: Vertical motion diagnosis
fig4 = figure('units','inch','position',[1,1,24,16], 'Visible', 'off');

subplot(2,3,1);
contourf(xx,zz,(squeeze(gpt_h(:,jj/2+1,:)))',-50:2:50,'linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca,'xtick',0:30:360);
set(gca,'ytick',0:1:10);
title('Cross section of geop. height at lat = 45')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,4);
ttt=squeeze(wfield(:,jj/2+1,:));
ttt(abs(ttt)<1.0e-8)=0; % set all noise data to zero
contourf(xx,zz,ttt','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca,'xtick',0:30:360);
set(gca,'ytick',0:1:10);
title('Cross section of vertical motions at lat = 45')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,2)
ttt=squeeze(temp(:,:,kk/2+1+10));
ttt(abs(ttt)<1.0e-8)=0;
contourf(xx,yy,ttt',-50:0.2:50,'linestyle', 'none');
hold on;
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1+10)))',0:2:50,'w');
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1+10)))',-50:2:-2,'--w');
colorbar;
set(gca,'clim',[-1 1]);
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Temp/height (shading/contour) at 300 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,5)
ttt=squeeze(temp(:,:,kk/2+1-10));
ttt(abs(ttt)<1.0e-8)=0; 
contourf(xx,yy,ttt',-50:0.2:50,'linestyle', 'none');
hold on;
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1-10)))',0:2:50,'w');
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1-10)))',-50:2:-2,'--w');
colorbar;
set(gca,'clim',[-1 1]);
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Temp/height (shading/contour) at 700 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,3)
contourf(xx,yy,(squeeze(wfield(:,:,kk/2+1+10)))','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vertical motion at 300 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,6)
contourf(xx,yy,(squeeze(wfield(:,:,kk/2+1-10)))','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vertical motion at 800 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

% Save the figure
ROW_VERT_MOTION = fullfile(row_plot_dir, 'row_vertical_motion_diagnosis.png');
saveas(fig4, ROW_VERT_MOTION);
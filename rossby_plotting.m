% rossby_plotting.m
% Load data, initialize variables, process variables, and perform plotting
load('data/Rossby_wave_2.mat')
addpath('functions');

global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx HH gg ii dx xx yy zz

% Initialize grid parameters
initialize_variables;

% Process variables
[eVec_amp, gpt_h, omega, phase_speed, growth_rate, eFolding, Amp, QV, XVy, XVx, XVz, temp, ug, vg, pvfield, gpt_h_hovmoler, pv, time] = process_variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1: Eigenvector amplitude
figure('units', 'inch', 'position', [4,2,16,12])
contourf(yy, zz, eVec_amp', 'linestyle', 'none');
xlabel('Latitude')
ylabel('Height (km)')
set(gca, 'fontsize', 22, 'color', 'w')
title =(['Zonal wave # = ', num2str(m0), ', eMode # =', num2str(7)...
    ', real eVal = ', num2str(growth_rate), ...
    ', imag eVal = ', num2str(omega)]);
colorbar

% Figure 2: Meridional cross-section of geopotential height
figure('units', 'inch', 'position', [4,2,16,12])
contourf(xx, zz, squeeze(gpt_h(:,jj/2+1,:))', 'LineStyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca, 'xtick', 0:30:360)
set(gca, 'ytick', 0:5:50)

% Figure 3: Hovmoller diagram of geopotential height
figure('units', 'inch', 'position', [4,2,16,12])
contourf(xx, time, gpt_h_hovmoler', -100:2:100, 'LineStyle', 'none');
set(gca,'xlim',[0,90])
colorbar;
xlabel('Longitude')
ylabel('Time (days)')
set(gca, 'xtick', 0:30:360)
set(gca, 'ytick', 0:5:50)
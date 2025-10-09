% initialize_variables.m
% Initialize constants and grid parameters
global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx HH gg ii dx xx yy zz

%% grid parameters
ii = 360; % longitude grid
dx = Lx/ii; % longitude grid spacing

xx = 0.0:360/ii:360; % longitude coordinates
yy = linspace(45-25, 45+25, jj+1); % latitude coordinates
zz = linspace(0.0,10,kk+1); % height grid (H=10km)
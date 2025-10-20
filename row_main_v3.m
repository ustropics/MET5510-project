%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: row_main.m

% Description: Script for setting up the Rossby wave model in the 
% quasi-geostrophic framework, initializing grid parameters, computing 
% the mean zonal wind (Ubar) and PV gradient (BPVy), constructing matrices 
% (B, C, D) for PV inversion and advection, and solving the eigenvalue 
% problem to determine wave stability, saving results to 
% 'data/rossby_wave_#.mat'.

% Functions used: 
% - l2jk: Converts linear index to 2D indices (j, k).
% - jk2l: Converts 2D indices (j, k) to linear index.
% - stream2pv: Computes potential vorticity from streamfunction.
% - stream2xPVadv: Calculates zonal advection of potential vorticity.
% - stream2yPVadv: Computes meridional advection of potential vorticity.
% - eig: Performs eigenvalue decomposition.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx
tic

addpath('functions'); % add functions folder
addpath('config'); % add config folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load constants from row_config.m and assign global variables
params = row_config();

% get directories
row_data_dir = params.row_data_dir;
row_data_filename = params.row_data_filename;
row_data = fullfile(row_data_dir, row_data_filename);

cplx = params.cplx; % imaginary unit for complex number operations
m0 = params.m0; % wave number
jj = params.jj; % number of latitude grid points
kk = params.kk; % number of height grid points
ll = params.ll; % total number of linear indices for interior points
Lx = params.Lx; % zonal domain length at 45° latitude (m)
f0 = params.f0; % coriolis parameter at 45° latitude (s^-1)
beta = params.beta; % beta parameter (meridional PV gradient, m^-1 s^-1)
gg = params.gg; % gravity (m s^-2)
Theta0 = params.Theta0; % reference potential temperature (K)
delta_Theta0 = params.delta_Theta0; % potential temperature difference (K)
HH = params.HH; % scale height (m)
dy = params.dy; % meridional grid spacing (in meters)
dz = params.dz; % vertical grid spacing (in meters)
NN2 = params.NN2; % brunt-vaisala frequency squared (s^-2)

%% Set Ubar with constant value
Ubar = zeros(jj+1, kk+1); % initialize mean zonal wind field
% Ubar(:,:) = Ubar(:,:) + params.Ubar_const; % set uniform zonal wind speed

% Ubar = hwe_ubar(params);

%% Initialize BPVy and set interior points to beta
% BPVy = zeros(jj+1, kk+1); % initialize beta-plane PV gradient
% for k = 2 : kk
%     for j = 2:jj
%         BPVy(j,k) = beta; % assign constant beta to the PV gradient
%     end
% end

%% Set up coordinates (already in your code, but included for context)
ZZ = 0.0:dz:HH;
Y0 = 6.37 * 1.0e6 * pi / 4;  % Central lat = 45°
Ly = 6.37 * 1.0e6 * 50 * pi / 180;
YY = Y0 - Ly/2 : dy : Y0 + Ly/2;

%% Parameters for Hoskins-West Ubar (from slide 9)
beta = 0;  % f-plane (b = 0)
Ubar0 = 0;  % Constant offset (can adjust if needed)
dTbar = 30;  % Align with delta_Theta0 = 30 K (or set to 40)
Ys = Y0;  % Central y for the jet
mu = 0.5;            % or 0 or 1 depending on test
Ly = 6.37e6 * 50 * pi/180;  % meridional domain (same as before)
Y0 = 6.37e6 * pi/4;  % center latitude (45°)
L_R = sqrt(NN2 * HH / f0^2);  % Rossby deformation radius
DeltaTbar = delta_Theta0;     % to match your variable names

% Preallocate
Ubar = zeros(jj+1, kk+1);

for j = 1:jj+1
    y = (j-1)*dy; % meridional coordinate
    for k = 1:kk+1
        z = (k-1)*dz; % vertical coordinate (0 to H)
        
        % Compute nondimensional z/H
        zH = z / HH;
        
        % Hoskins-West mean flow
        term1 = zH - (mu/2) * ( zH + ...
            (sinh((2*pi*L_R*z)/(Ly*HH)) / sinh((2*pi*L_R)/Ly)) * ...
            cos(2*pi*(y - Y0)/Ly) );
        
        Ubar(j,k) = (gg * HH * DeltaTbar) / (f0 * Theta0 * Ly) * term1;
    end
end


BPVy = zeros(jj+1, kk+1);
for j = 2:jj
    for k = 2:kk
        BPVy(j,k) = beta - (Ubar(j+1,k)-2*Ubar(j,k) + Ubar(j-1,k)) / dy^2 ...
            -(f0*f0/NN2) * (Ubar(j,k+1) - 2*Ubar(j,k) + Ubar(j,k-1))/dz^2;
    end

    k=1;
    BPVy(j,k)=-(Ubar(j,k+1) - Ubar(j,k))/dz;
    k = kk+1;
    BPVy(j,k)=-(Ubar(j,k) - Ubar(j,k-1))/dz;
end

%% This is just for testing
[j,k] = l2jk(98); % convert linear index 98 to 2D indices (j,k)
y = jk2l(j,k); % convert back to linear index for testing

%% Initialize vectors for PV and its advections
XV = zeros(ll,1); % streamfunction vector
QV = zeros(ll,1); % PV vector
xQVadv = zeros(ll,1); % zonal advection of PV
yQVadv = zeros(ll,1); % meridional advection of PV

%% Initialize matrices
B = zeros(ll,ll); % matrix for PV inversion
C = zeros(ll, ll); % matrix for zonal PV advection
D = zeros(ll,ll); % matrix for meridional PV advection (beta term = constant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP & EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main loop to construct matrices B, C, and D
for l0 = 1:ll
    XV = zeros(ll,1); % reset streamfunction vector
    XV(l0) = 1;  % set to 1 for matrix construction

    % B MATRIX
    QV = stream2pv(XV);
    B(:, l0) = QV(:);

    % C MATRIX
    xQVadv = stream2xPVadv(QV);
    C(:, l0) = xQVadv(:);

    % D MATRIX
    yQVadv = stream2yPVadv(XV);
    D(:,l0) = yQVadv(:);
end

%% Start eigenvector and eigenvalue part
% Linearized PV equation in matrix form for eigenvalue problem
A = B^(-1) * (C+D);

[eigVec, eigValm] = eig(A); % eigVec = eigenvectors
eigVal = diag(eigValm); % eigValm = eigenvalue matrix

% Sort by descending real part
[test, sortdx] = sort(real(eigVal), 'descend');

eigVal2 = eigVal(sortdx);
eigVec2 = eigVec(:, sortdx);

% Sort by ascending real part of imaginary eigenvalue
[test, sortdx] = sort(real(cplx*eigVal), 'ascend');
eigVal3 = eigVal(sortdx);
eigVec3 = eigVec(:,sortdx);

toc

%% Create data folder if it doesn't exist and save data file
if ~exist(row_data_dir, 'dir')
    mkdir(row_data_dir);
end

% save it
save(row_data, 'Ubar', 'BPVy', 'eigVec', 'eigVal', 'eigVec2', 'eigVal2', ...
    'eigVec3', 'eigVal3', 'jj', 'kk', 'll', 'NN2', 'f0', 'dy', 'dz', 'm0', ...
    'Lx', 'beta', 'cplx', 'gg', 'Theta0', 'HH');
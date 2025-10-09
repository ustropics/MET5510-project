global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These can be found in Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% beginning on slides 5 - 6

%% Initial variables for complex and wave number
cplx = sqrt(-1); % imaginary unit for complex number operations
m0 = 2; % this is wave number (0-3 is longwave/planetary, 3+ is shortwave)

%% grid point parameterization 
jj = 50; % number of latitude grid points
kk = 50; % number of height grid points
ll = (jj-1)*(kk+1); % total number of linear indices for interior points
Lx = 2 * pi * 6.37 * 1.0e6 * cos(pi/4); % zonal domain length at 45° latitude (m)

%% Coriolis, beta, and gravity calculations/constants
f0 = 2 * (2*pi/86400) * sin(pi/4); % coriolis parameter at 45° latitude (s^-1)
beta = 2 * (2 * pi/86400) * cos(pi/4)/(6.37*1.e6); % beta parameter (meridional PV gradient, m^-1 s^-1)
gg = 9.81; % gravity (m s^-2)

%% Reference potential temperatures, height, and additional domain gridding
Theta0 = 300; % reference potential temperature (K)
delta_Theta0 = 30; % potential temperature difference across dmain (K)
HH = 10000; % scale height (m)

Ly = 6.37 * 1.0e6*50 * pi/180; % meridional domain length (50° latitude range, m)
dy = Ly/jj; % meridional grid spacing (in meters)
dz = HH/kk; % vertical grid spacing (in meters)

%% Set brunt-vaisala frequency, ubar, and beta-plane (constant)
NN2 = gg * delta_Theta0/(HH*Theta0); % brunt-vaisala frequency squared (s^-2)
Ubar = zeros(jj+1,kk+1); % initialize mean zonal wind field (constant 10 m s^-1)
Ubar(:,:) = Ubar(:,:) + 10; % set uniform zonal wind speed
BPVy = zeros(jj+1,kk+1); % initialize beta-plane PV gradient

%% This is just for testing
[j,k] = l2jk(98); % convert linear index 98 to 2D indices (j,k)
y = jk2l(j,k); % convert back to linear index for testing

%% set BPVy to beta for interior points
for k = 2 : kk
    for j = 2:jj
        BPVy(j,k) = beta; % assign constant beta to the PV gradient
    end
end

%% Initialize vectors for PV and its advections
XV = zeros(ll,1); % streamfunction vector
QV = zeros(ll,1); % PV vector
xQVadv = zeros(ll,1); % zonal advection of PV
yQVadv = zeros(ll,1); % meridional advection of PV

%% Initialize matrices
B = zeros(ll,ll); % matrix for PV inversion
C = zeros(ll, ll); % matrix for zonal PV advection
D = zeros(ll,ll); %  matrix for meridional PV advection (beta term = constant)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP & EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main loop to construct matrices B, C, and D
% This is from page 16 from Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% (Four vectors involved in the model)
for l0 = 1:ll
    XV = zeros(ll,1); % reset streamfunction vector
    XV(l0) = 1;  % set 10m to level 1

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

% test values in descending order
[test, sortdx] = sort(real(eigVal), 'descend');


% this is a solution to the real part
eigVal2 = eigVal(sortdx);
eigVec2 = eigVec(:, sortdx);

% test values in ascending order
[test, sortdx] = sort(real(cplx*eigVal), 'ascend');

% this is a solution for the complex/imaginary part
eigVal3 = eigVal(sortdx);
eigVec3 = eigVec(:,sortdx);

toc

save('Rossby_wave_2.mat') % save data file for use later on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% jk2l function: identifies 2d (j,k) indices to the grid point
function l = jk2l(j,k)
global jj

l = j - 1 + (k-1) * (jj-1);

end

%% l2jk function: reverse of jk2l, identifies grid point to 2d (j,k) indices
function [j,k] = l2jk(l)
global jj

k = floor((l-1)/(jj-1))+1; % vertical index
j = l+1 - (k-1)*(jj-1); % meridional index
end

%% QV function: This inverts the PV-streamfunction relation using finite differences
% This can be found in Dr. Cai's 3D_spectral_linear_QG_model_numnerics.pdf
% for combining equations in interior and at the top/boundary (page 8)
% this also 

function QV = stream2pv(XV)
global jj kk ll NN2 m0 f0 dy dz Lx

QV = zeros(ll,1);

% bottom boundary (k = 1)
k = 1; 
for j = 2:jj
    l = jk2l(j,k); % from j value we can get eigen value
    lup = jk2l(j, k+1); % get one level from above
    ldn = l; % we've reached bottom boundary

    QV(l) = ( XV(lup) - XV(ldn) )/dz;
end

% top boundary (k=kk+1)
k = kk + 1;
for j = 2:jj
    l = jk2l(j, k);
    ldn = jk2l(j, k-1);
    lup = l; % we've reached top boundary

    QV(l) = ( XV(lup) - XV(ldn) ) / dz;
end

% interior points (quasi-geostrophic PV)
for k = 2:kk

    for j = 2:jj

        l = jk2l(j, k); % current index
        lup = jk2l(j, k+1); % above
        ldn = jk2l(j, k-1); % below
        lnh = jk2l(j+1,k); % north (+j)
        lsh = jk2l(j-1,k); % south (-j)

        % meridional boundary condition
        if (j == 2)
            XVsh=0; % southern boundary
        else
            XVsh = XV(lsh);
        end

        if (j == jj)
            XVnh = 0; % northern boundary

        else
            XVnh = XV(lnh);

        end

        % PV = - (k^2 + laplacian_psi) + (f0^2 / N^2) * d^2 psi / dz^2
        % Where k = 2*pi*m0/Lx  (zonal wavenumber)
        QV(l) = -((2*pi*m0/Lx)^2)*XV(l) ...
            + ( XVsh - 2*XV(l) + XVnh) / dy/dy ... 
            + (f0/dz)^2 * (XV(lup) - 2*XV(l) + XV(ldn)) / NN2;

    end

end
end

%% Function xQVadv: compute zonal advection of PV
% -Ubar * ik * Q, where k is wavenumber
function xQVadv = stream2xPVadv(QV)
global jj kk ll Lx Ubar cplx m0

xQVadv = zeros(ll,1);

for k = 1:kk + 1

    for j = 2:jj

        l = jk2l(j,k);

        xQVadv(l) = -Ubar(j,k) * cplx * (2*pi*m0/Lx) * QV(l);

    end
end
end

%% Function yQVadv: compute meridional advection of PV
% -Ubar * ik * Q, where k is wavenumber
function yQVadv = stream2yPVadv(XV)
global jj kk ll BPVy Lx cplx m0

yQVadv = zeros(ll,1);

for k = 1:kk + 1

    for j = 2:jj

        l = jk2l(j,k);

        yQVadv(l) = -BPVy(j,k) * cplx * (2*pi*m0/Lx) * XV(l);
    end
end
end
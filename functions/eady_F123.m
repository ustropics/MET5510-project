%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: eady_F123.m

% DESCRIPTION: Computes the forcing terms F1, F2, and F3 for vertical motion 
% calculations in the Eady wave model within the quasi-geostrophic framework. 
% These terms represent contributions to the vertical velocity equation, 
% incorporating the effects of zonal wind shear, potential vorticity gradients, 
% and beta-plane dynamics. The function is designed to be used with the 
% eady_plot.m script for generating vertical motion diagnostics.

% INPUT:
% - jj: Number of latitude grid points (integer)
% - kk: Number of height grid points (integer)
% - dy: Grid spacing in the meridional direction (m)
% - dz: Grid spacing in the vertical direction (m)
% - m0: Zonal wavenumber (integer)
% - Lx: Zonal domain length (m)
% - cplx: Imaginary unit for complex number operations (sqrt(-1))
% - beta: Planetary vorticity gradient (s^-1 m^-1)
% - XV: Streamfunction vector (size ll x 1)
% - Ubar: Mean zonal wind field (jj+1 x kk+1 array, m/s)

% OUTPUT:
% - F1: Forcing term related to vorticity advection (LW x 1 vector)
% - F2: Forcing term related to meridional PV gradient (LW x 1 vector)
% - F3: Forcing term related to beta effect (LW x 1 vector)

% MATH/FUNCTIONS: 
% - F1: Represents vorticity advection term, computed as:
%       F1 = 2 * (2πm₀/Lx) * i * (∂U/∂z) * [-(2πm₀/Lx)²ψ + ∂²ψ/∂y²] / (2dz)
% - F2: Represents meridional PV gradient term, computed as:
%       F2 = -2 * (2πm₀/Lx) * i * (∂²U/∂y²) * (∂ψ/∂z) / (2dz * dy²)
% - F3: Represents beta effect term, computed as:
%       F3 = (2πm₀/Lx) * i * β * (∂ψ/∂z) / (2dz)

% - where:
%   - m₀: Zonal wavenumber
%   - Lx: Zonal domain length
%   - i: Imaginary unit (cplx)
%   - β: Planetary vorticity gradient
%   - U: Mean zonal wind (Ubar)
%   - ψ: Streamfunction (XV)
%   - dy, dz: Grid spacings
% - Uses jk2lw to map 2D indices (j,k) to linear weighted indices for vertical motion.

% DEPENDENCIES:
% - jk2lw: Converts 2D indices (j, k) to a weighted linear index.
% - jk2l: Converts 2D indices (j, k) to a linear index.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F1, F2, F3] = eady_F123(params, XV, Ubar)

    %% Extract parameters
    jj = params.jj;
    kk = params.kk;
    dy = params.dy;
    dz = params.dz;
    m0 = params.m0;
    Lx = params.Lx;
    cplx = params.cplx;
    beta = params.beta;
    
    %% Initialize output arrays
    LW = (jj+1)*(kk-1); % Size for vertical motion matrix
    F1 = zeros(LW,1);
    F2 = zeros(LW,1);
    F3 = zeros(LW,1);

    %% Compute F2 and F3 for interior latitude points (j = 2 to jj)
    for j = 2:jj
        for k = 2:kk
            l = jk2lw(j,k);
            F3(l) = (2*pi*m0/Lx)*cplx*beta*(XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz);
            F2(l) = -2*(2*pi*m0/Lx)*cplx*(Ubar(j+1,k)-2*Ubar(j,k)+Ubar(j-1,k)) ...
                    *(XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz*dy*dy);
        end
    end

    %% Compute F1 for boundary latitude (j = 1)
    j = 1;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
            *(XV(jk2l(3,k))-2*XV(jk2l(2,k)))/(2*dz*dy*dy);
    end

    %% Compute F1 for latitude j = 2
    j = 2;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
            *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
            (XV(jk2l(3,k))-2*XV(jk2l(2,k)))/dy/dy)/(2*dz);
    end

    %% Compute F1 for latitude j = jj
    j = jj;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
            *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
             (XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/dy/dy)/(2*dz);
    end

    %% Compute F1 for latitude j = jj + 1
    j = jj + 1;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
              *(XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/(2*dz*dy*dy);
    end

    %% Compute F1 for interior latitudes (j = 3 to jj-1)
    for j = 3:jj-1
        for k = 2:kk
            l = jk2lw(j,k);
            F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
                *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
                 (XV(jk2l(j+1,k))-2*XV(jk2l(j,k))+XV(jk2l(j-1,k)))/dy/dy)/(2*dz);
        end
    end
end
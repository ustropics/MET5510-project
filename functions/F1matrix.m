%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: F1matrix.m
%
% DESCRIPTION: Computes the vertical shear forcing term (F1) in the QG omega 
% equation, representing geostrophic deformation acting on background vertical 
% shear. Uses specialized finite-difference stencils at meridional boundaries 
% and interior points to maintain accuracy and boundary condition consistency.
%
% INPUT:
% - XV: 1D streamfunction eigenvector (m²/s), length ll
% - Ubar: 2D background zonal wind (m/s), size (jj+1, kk+1)
%
% OUTPUT:
% - F1: 1D forcing vector (s⁻²), length LW
%
% MATH/FUNCTIONS: 
% - F1 = 2 i k (∂U/∂z) [ -k²ψ + ∂²ψ/∂y² ] / (2 dz)
%
% VARIABLES:
% - k = 2π m0 / Lx: zonal wavenumber
% - ∂U/∂z: centered in z: (U(k+1) - U(k-1))/(2 dz)
% - At j=1, jj+1 (walls): uses one-sided y-stencil on ψ
% - At j=2, jj: uses ghost-node-aware centered stencil
% - At interior 3 ≤ j ≤ jj-1: full centered second y-derivative
% - m0, Lx, dy, dz, cplx: model parameters
% - jk2l, jk2lw: index mapping functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F1 = F1matrix(XV, Ubar)
    global jj kk cplx m0 Lx dy dz

    LW = (jj+1)*(kk-1);
    F1 = zeros(LW, 1);

    % Precompute common factor: 2 i k = 2 * cplx * (2π m0 / Lx)
    ik2 = 2 * cplx * (2*pi*m0/Lx);

    %% j = 1 (southern boundary, j=1 is ghost-adjacent)
    j = 1;
    for k = 2:kk
        l = jk2lw(j, k);
        dUdz = (Ubar(j,k+1) - Ubar(j,k-1)) / (2*dz);  % centered in z
        d2psidy2 = (XV(jk2l(3,k)) - 2*XV(jk2l(2,k))) / (dy*dy);  % one-sided
        F1(l) = ik2 * dUdz * d2psidy2 / (2*dz);
    end

    %% j = 2 (first interior, uses ghost at j=1=0)
    j = 2;
    for k = 2:kk
        l = jk2lw(j, k);
        dUdz = (Ubar(j,k+1) - Ubar(j,k-1)) / (2*dz);
        k2psi = -(2*pi*m0/Lx)^2 * XV(jk2l(j,k));
        d2psidy2 = (XV(jk2l(3,k)) - 2*XV(jk2l(2,k))) / (dy*dy);
        F1(l) = ik2 * dUdz * (k2psi + d2psidy2) / (2*dz);
    end

    %% j = jj (last interior, uses ghost at j=jj+1=0)
    j = jj;
    for k = 2:kk
        l = jk2lw(j, k);
        dUdz = (Ubar(j,k+1) - Ubar(j,k-1)) / (2*dz);
        k2psi = -(2*pi*m0/Lx)^2 * XV(jk2l(j,k));
        d2psidy2 = (XV(jk2l(jj-1,k)) - 2*XV(jk2l(jj,k))) / (dy*dy);
        F1(l) = ik2 * dUdz * (k2psi + d2psidy2) / (2*dz);
    end

    %% j = jj+1 (northern boundary)
    j = jj+1;
    for k = 2:kk
        l = jk2lw(j, k);
        dUdz = (Ubar(j,k+1) - Ubar(j,k-1)) / (2*dz);
        d2psidy2 = (XV(jk2l(jj-1,k)) - 2*XV(jk2l(jj,k))) / (dy*dy);
        F1(l) = ik2 * dUdz * d2psidy2 / (2*dz);
    end

    %% interior latitudes 3 ≤ j ≤ jj-1
    for j = 3:jj-1
        for k = 2:kk
            l = jk2lw(j, k);
            dUdz = (Ubar(j,k+1) - Ubar(j,k-1)) / (2*dz);
            k2psi = -(2*pi*m0/Lx)^2 * XV(jk2l(j,k));
            d2psidy2 = (XV(jk2l(j+1,k)) - 2*XV(jk2l(j,k)) + XV(jk2l(j-1,k))) / (dy*dy);
            F1(l) = ik2 * dUdz * (k2psi + d2psidy2) / (2*dz);
        end
    end
end
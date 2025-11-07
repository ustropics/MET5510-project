%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: F2matrix.m
%
% DESCRIPTION: Computes the horizontal shear forcing term (F2) in the QG omega 
% equation, representing the contribution from differential zonal advection of 
% background vorticity by the perturbation. Active only at interior meridional 
% points due to boundary conditions.
%
% INPUT:
% - XV: 1D array representing the streamfunction eigenvector (m²/s), length ll
% - Ubar: 2D background zonal wind (m/s), size (jj+1, kk+1)
%
% OUTPUT:
% - F2: 1D forcing vector (s⁻²), length LW (interior points only)
%
% MATH/FUNCTIONS: 
% - F2 = -2 i k (∂U/∂y) (∂²ψ/∂y ∂z) / (2 dz dy²)
%
% VARIABLES:
% - k = 2π m0 / Lx: zonal wavenumber (m⁻¹)
% - ∂U/∂y: centered second difference in y: (U(j+1)-2U(j)+U(j-1))/dy²
% - ∂²ψ/∂y ∂z: mixed derivative using centered differences
% - m0: Zonal wavenumber (dimensionless)
% - Lx: Zonal domain length (m)
% - cplx: Complex unit i = √(-1)
% - dy, dz: Grid spacing in y, z (m)
% - jk2l: Maps (j,k) → linear index in XV (ll-sized)
% - jk2lw: Maps (j,k) → linear index in F2/w_vec (LW-sized)
% - Applied only for 2 ≤ j ≤ jj, 2 ≤ k ≤ kk (interior)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F2 = F2matrix(XV, Ubar)
    global jj kk cplx m0 Lx dy dz
    
    % Preallocate F2 vector for interior points only (LW = (jj+1)*(kk-1))
    F2 = zeros((jj+1)*(kk-1), 1);

    % Loop over interior meridional (j) and vertical (k) indices
    for j = 2:jj
        for k = 2:kk
            % Linear index in output vector F2 (LW-sized)
            l = jk2lw(j, k);
            
            % Zonal wavenumber factor: i * k = i * (2π m0 / Lx)
            ik = cplx * (2*pi*m0/Lx);
            
            % Second derivative of Ubar in y: ∂²U/∂y² ≈ (U(j+1) - 2U(j) + U(j-1))/dy²
            d2Udy2 = (Ubar(j+1,k) - 2*Ubar(j,k) + Ubar(j-1,k)) / (dy*dy);
            
            % Mixed derivative ∂²ψ/∂y∂z ≈ [ψ(j,k+1) - ψ(j,k-1)] / (2 dz)
            d2psidydz = (XV(jk2l(j,k+1)) - XV(jk2l(j,k-1))) / (2*dz);
            
            % Full F2 term: -2 i k (∂²U/∂y²) (∂²ψ/∂y∂z) / (2 dz)
            % Note: the extra /dy² comes from ∂²U/∂y²
            F2(l) = -2 * ik * d2Udy2 * d2psidydz / (2*dz);
        end
    end
end
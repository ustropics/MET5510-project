%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: F3matrix.m
%
% DESCRIPTION: Computes the planetary vorticity advection term (F3) in the QG 
% omega equation, representing the beta effect on vertical motion. Only active 
% at interior meridional grid points due to boundary treatment.
%
% INPUT:
% - XV: 1D streamfunction eigenvector (m²/s), length ll
%
% OUTPUT:
% - F3: 1D forcing vector (s⁻²), length LW
%
% MATH/FUNCTIONS: 
% - F3 = i k β (∂ψ/∂z) / (2 dz)
%
% VARIABLES:
% - k = 2π m0 / Lx: zonal wavenumber
% - β: Planetary vorticity gradient (m⁻¹ s⁻¹)
% - ∂ψ/∂z: centered vertical derivative: (ψ(k+1) - ψ(k-1))/(2 dz)
% - m0, Lx, dz, cplx: model parameters
% - jk2l, jk2lw: index mapping
% - Applied only for 2 ≤ j ≤ jj, 2 ≤ k ≤ kk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F3 = F3matrix(XV)
    global jj kk cplx m0 Lx beta dz

    % Preallocate F3 for interior points
    F3 = zeros((jj+1)*(kk-1), 1);

    % Precompute i k = cplx * (2π m0 / Lx)
    ik = cplx * (2*pi*m0/Lx);

    % Loop over interior (j,k) points only
    for j = 2:jj
        for k = 2:kk
            % Linear index in F3 vector
            l = jk2lw(j, k);
            
            % Vertical derivative of streamfunction: ∂ψ/∂z
            dpsi_dz = (XV(jk2l(j,k+1)) - XV(jk2l(j,k-1))) / (2*dz);
            
            % Full F3 term: i k β (∂ψ/∂z) / (2 dz)
            F3(l) = ik * beta * dpsi_dz / (2*dz);
        end
    end
end
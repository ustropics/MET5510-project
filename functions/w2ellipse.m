%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: w2ellipse.m

% DESCRIPTION: Applies the elliptic diagnostic operator to the vertical velocity 
% field to recover the forcing balance in the omega equation of the 
% quasi-geostrophic (QG) model. This function computes the elliptical 
% representation used to validate the inversion of the QG omega equation.

% INPUT:
% - w: 1D array representing the interior vertical velocity (m/s), with length LW

% OUTPUT:
% - EW: 1D array representing the elliptic operator applied to w (s⁻²), 
%        with length LW

% MATH/FUNCTIONS: 
% - EW = -k²w + ∂²w/∂y² + (f₀²/N²) ∂²w/∂z²

% VARIABLES:
% - k = 2π m0 / Lx: zonal wavenumber (m⁻¹)
% - ∂²w/∂y²: second meridional derivative via centered finite differences
% - ∂²w/∂z²: second vertical derivative via centered/one-sided differences
% - m0: Zonal wavenumber (dimensionless)
% - Lx: Zonal domain length (m)
% - f0: Coriolis parameter (s⁻¹)
% - NN2: Brunt-Väisälä frequency squared, N² (s⁻²)
% - dy: Meridional grid spacing (m)
% - dz: Vertical grid spacing (m)
% - jk2lw: Function mapping (j,k) → linear index in w_vec (LW-sized)
% - The operator is discretized using second-order finite differences, 
%   with one-sided stencils at near-boundary points and Dirichlet w=0 at walls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EW = w2ellipse(w)
    global jj kk LW NN2 m0 f0 dy dz Lx

    EW = zeros(LW,1);

    %% Southern wall (j = 1)
    j = 1; % initialize j for northern wall loop
    for k = 2:kk
        l   = jk2lw(j,k);
        ln3 = jk2lw(3,k);
        ln2 = jk2lw(2,k);

        % Handle vertical boundaries
        if k == 2
            wdn = 0; wup = w(jk2lw(j,k+1)); % bottom boundary
        elseif k == kk
            wdn = w(jk2lw(j,k-1)); wup = 0; % top boundary

        % Interior points
        else
            wdn = w(jk2lw(j,k-1)); wup = w(jk2lw(j,k+1)); % interior
        end

        % Compute elliptic operator
        EW(l) = -(2*pi*m0/Lx)^2*w(l) ...
                + (w(ln3)-2*w(ln2)+w(l))/dy^2 ...
                + (f0/dz)^2*(wup-2*w(l)+wdn)/NN2;
    end

    %% Northern wall (j = jj+1)
    j = jj+1; % initialize j for northern wall loop
    for k = 2:kk
        l   = jk2lw(j,k);
        ls3 = jk2lw(jj-1,k);
        ls2 = jk2lw(jj,k);

        % Handle vertical boundaries
        if k == 2
            wdn = 0; wup = w(jk2lw(j,k+1)); % bottom boundary
        elseif k == kk
            wdn = w(jk2lw(j,k-1)); wup = 0; % top boundary
        else
            wdn = w(jk2lw(j,k-1)); wup = w(jk2lw(j,k+1)); % interior
        end

        % Compute elliptic operator
        EW(l) = -(2*pi*m0/Lx)^2*w(l) ...
                + (w(l)-2*w(ls2)+w(ls3))/dy^2 ...
                + (f0/dz)^2*(wup-2*w(l)+wdn)/NN2;
    end

    %% Interior points
    for j = 2:jj
        for k = 2:kk % loop over interior grid points
            l   = jk2lw(j,k);
            ln1 = jk2lw(j+1,k);
            ls1 = jk2lw(j-1,k);

            % Handle vertical boundaries
            if k == 2
                wdn = 0; wup = w(jk2lw(j,k+1)); % bottom boundary
            elseif k == kk
                wdn = w(jk2lw(j,k-1)); wup = 0; % top boundary
            else
                wdn = w(jk2lw(j,k-1)); wup = w(jk2lw(j,k+1)); % interior
            end

            % Compute elliptic operator
            EW(l) = -(2*pi*m0/Lx)^2*w(l) ...
                    + (w(ln1)-2*w(l)+w(ls1))/dy^2 ...
                    + (f0/dz)^2*(wup-2*w(l)+wdn)/NN2;
        end
    end
    
end
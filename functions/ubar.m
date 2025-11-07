%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: ubar.m

% DESCRIPTION: Constructs the background zonal mean wind field Ubar(y,z) for the 
% Hoskins-West modified Eady (HWME) model in the quasi-geostrophic framework. 
% The profile includes a linear shear in height, a curvature term controlled by mu, 
% and a meridionally varying vertical structure via a hyperbolic sine function, 
% designed to produce a baroclinic jet with realistic thermal wind balance.

% INPUT:
% - jj: Number of meridional grid points (interior-adjacent)
% - kk: Number of vertical grid points (including boundaries)
% - gg: Gravitational acceleration (m s⁻²)
% - f0: Coriolis parameter at reference latitude (s⁻¹)
% - Theta0: Reference potential temperature (K)
% - dTbar: Horizontal potential temperature gradient across domain (K)
% - HH: Scale height (m)
% - Ly: Meridional domain length (m)
% - ZZ: 1D array of vertical coordinates (m), length kk+1
% - YY: 1D array of meridional coordinates (m), length jj+1
% - mu: Jet curvature parameter (dimensionless)
% - Lr: Internal Rossby deformation radius, N H / f0 (m)

% OUTPUT:
% - Ubar: 2D array of mean zonal wind (m/s), size (jj+1, kk+1)

% MATH/FUNCTIONS: 
%   U(y,z) = prefac * [ (z/H) - (mu/2)(z/H)² + (sinh(λ z/H)/sinh(λ)) * cos(2π (y-y_s)/Ly) ]
%   where prefac = (g / (f0 Theta0)) * (H dTbar / Ly),   lambda = 2π Lr / Ly

% VARIABLES:
% - prefac: Thermal wind scaling factor (m/s)
% - lambda = 2π Lr / Ly: nondimensional vertical structure parameter
% - zH = z/H: normalized height
% - y_s = y - y_min: distance from southern boundary (m)
% - sinh_term = sinh(lambda zH) / sinh(lambda): vertical structure function
% - cos_term = cos(2π (y-y_s)/Ly): meridional modulation
% - The profile satisfies thermal wind balance with a uniform horizontal 
%   temperature gradient dTbar and includes a Gaussian-like jet via mu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ubar = ubar(jj, kk, gg, f0, Theta0, dTbar, HH, Ly, ZZ, YY, mu, Lr)
    
    Ubar = zeros(jj + 1, kk + 1);
    prefac = (gg / (f0 * Theta0)) * (HH * dTbar / Ly); % overall prefactor
    param = 2 * pi * Lr / Ly; % parameter for vertical structure
    
    %% Compute Ubar(y,z) over the grid
    for j = 1:(jj + 1) % loop over meridional grid points
        y_s = (YY(j)-YY(1)); % meridional distance from southern boundary

        % Loop over vertical grid points
        for k = 1:(kk + 1)
            zH = ZZ(k) / HH; % normalized height
            sinh_term = sinh(param * zH) / sinh(param); % vertical structure term
            cos_term = cos(2 * pi * (y_s) / Ly); % meridional modulation term

            % Compute Ubar at (j, k)
            Ubar(j, k) = prefac * (zH - 0.5 * mu * (zH + sinh_term * cos_term)); 
        end
   
    end

end
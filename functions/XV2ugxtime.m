%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XV2ugxtime.m

% DESCRIPTION: Computes the time evolution of the zonal wind (ug) along longitude
% for a Hovmoller diagram at a specific latitude and height.

% INPUT:
% - XVy: Meridional derivative of streamfunction vector (ll x 1)
% - ii: Number of longitude grid points
% - dx: Longitudinal grid spacing (m)
% - omega: Imaginary part of eigenvalue (for time evolution)
% - ylat: Latitude index for Hovmoller
% - zlev: Vertical level index for Hovmoller

% OUTPUT:
% - ugxtime: Hovmoller data for zonal wind (ii+1 x 51 array, m/s)

% DEPENDENCIES:
% - XV2field: Computes 3D field from streamfunction vector
% - jk2l: Converts 2D indices (j, k) to linear index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ugxtime = XV2ugxtime(XVy, ii, dx, omega, ylat, zlev)
global cplx m0 Lx

ugxtime = zeros(ii+1, 51);

l = jk2l(ylat, zlev);
for day = 1:51
    sec = (day-1)*86400; % time in seconds
    for n = 1:ii+1 % loop over longitude grid points
        xlon = (n-1)*dx;

        % Compute ug at (n, day)
        ugxtime(n, day) = real( - XVy(l) * exp(cplx * (2 * pi * m0 * xlon / Lx + omega * sec)));
    end
end
end
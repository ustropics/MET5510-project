%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2ugxtime.m

% Description: Computes the time evolution of the zonal wind (ug) along longitude
% for a Hovmoller diagram at a specific latitude and height.

% Input:
% - XVy: Meridional derivative of streamfunction vector (ll x 1)
% - ii: Number of longitude grid points
% - dx: Longitudinal grid spacing (m)
% - omega: Imaginary part of eigenvalue (for time evolution)
% - ylat: Latitude index for Hovmoller
% - zlev: Vertical level index for Hovmoller

% Output:
% - ugxtime: Hovmoller data for zonal wind (ii+1 x 51 array, m/s)

% Dependencies:
% - XV2field: Computes 3D field from streamfunction vector
% - jk2l: Converts 2D indices (j, k) to linear index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ugxtime = XV2ugxtime(XVy, ii, dx, omega, ylat, zlev)
global cplx m0 Lx

ugxtime = zeros(ii+1, 51);

% Compute the full ug field at the specified latitude and height
% ug_field = -XV2field(XVy, ii, dx); % 3D field (ii+1 x jj+1 x kk+1)
% ug_slice = squeeze(ug_field(:, ylat, zlev)); % Extract 1D longitude slice

for day = 1:51
    sec = (day-1)*86400;
    for n = 1:ii+1
        xlon = (n-1)*dx;
        ugxtime(n, day) = real(XVy(ylat,zlev) * exp(cplx * (2 * pi * m0 * xlon / Lx + omega * sec)));
    end
end
end
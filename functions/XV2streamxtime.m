%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XV2streamxtime.m

% DESCRIPTION: Converts the 1D eigenvector (XV) into a streamfunction field 
% varying with longitude and time, tracking wave propagation in 
% the quasi-geostrophic model.

% INPUT:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction
% - omega: Growth rate from eigenvalues
% - ylat: Latitude index
% - zlev: Height index

% OUTPUT:
% - xtime: 2D array representing the streamfunction field over longitude and time

% MATH/FUNCTIONS: 
% - ψ(x,t) = XV * exp(i(kx - ωt))

% VARIABLES:
% - k is wavenumber
% - ω is growth rate from eigenvalues
% - x is longitude
% - t is time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xtime = XV2streamxtime(XV,ii,dx,omega,ylat,zlev)

    global cplx m0 Lx
    
    xtime = zeros(ii+1, 51); % initialize output array
    
    %% Get linear index for specified latitude and height
    l = jk2l(ylat, zlev);
    for day = 1: 51 % loop over days
        sec = (day-1)*86400; % time in seconds
        for n = 1:ii+1 % loop over longitude grid points
            xlon=(n-1)*dx; % longitudinal position

            % Compute streamfunction at (n, day)
            xtime(n,day) = real(XV(l)*exp(cplx*(2*pi*m0*xlon/Lx+omega*sec)));
        end
    end

end
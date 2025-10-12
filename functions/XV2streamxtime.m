%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2streamxtime.m

% Description: Converts the 1D eigenvector (XV) into a streamfunction field 
% varying with longitude and time, tracking wave propagation in 
% the quasi-geostrophic model.

% Input:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction
% - omega: Growth rate from eigenvalues
% - ylat: Latitude index
% - zlev: Height index

% Output:
% - xtime: 2D array representing the streamfunction field over longitude and time

% Math/functions: ψ(x,t) = XV * exp(i(kx - ωt)), where
% k is wavenumber
% ω is growth rate from eigenvalues
% x is longitude
% t is time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xtime = XV2streamxtime(XV,ii,dx,omega,ylat,zlev)
    global cplx m0 Lx
    
    xtime = zeros(ii+1, 51);
    
    l = jk2l(ylat, zlev);
    for day = 1: 51
        sec = (day-1)*86400;
        for n = 1:ii+1
            xlon=(n-1)*dx;
            xtime(n,day) = real(XV(l)*exp(cplx*(2*pi*m0*xlon/Lx+omega*sec)));
        end
    end
    end
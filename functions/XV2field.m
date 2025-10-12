%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2field.m

% Description: Transforms the 1D eigenvector (XV) into a 3D field 
% (e.g., geopotential height) over longitude, latitude, and height grids 
% using spectral or finite difference methods in the quasi-geostrophic model.

% Input:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction

% Output:
% - field: 3D array representing the geopotential height field

% Math/functions: Field = XV * exp(i * k * x), where 
% XV is the eigenvector
% k is the zonal wavenumber
% x is the longitude coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field= XV2field(XV,ii,dx)
    global jj kk ll cplx m0 Lx
    
    field=zeros(ii+1,jj+1,kk+1);
    
    for l = 1:ll
        [j,k]=l2jk(l);
        for i = 1:ii+1
            xlon=(i-1)*dx;
            % compute real part of complex exponential 
            field(i,j,k)=real(XV(l)*exp(cplx*2*pi*m0*xlon/Lx));
        end
    end
    end
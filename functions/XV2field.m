%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XV2field.m

% DESCRIPTION: Transforms the 1D eigenvector (XV) into a 3D field 
% (e.g., geopotential height) over longitude, latitude, and height grids 
% using spectral or finite difference methods in the quasi-geostrophic model.

% INPUT:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction

% OUTPUT:
% - field: 3D array representing the geopotential height field

% MATH/FUNCTIONS: Field = XV * exp(i * k * x)

% - VARIABLES:
%   - XV is the eigenvector
%   - k is the zonal wavenumber
%   - x is the longitude coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field= XV2field(XV,ii,dx)

    global jj kk ll cplx m0 Lx
    
    field=zeros(ii+1,jj+1,kk+1); % initialize output array
    
    for l = 1:ll
        [j,k]=l2jk(l); % convert linear index to (j,k)
        for i = 1:ii+1 % loop over longitude grid points
            xlon=(i-1)*dx; % longitudinal position

            % compute real part of complex exponential 
            field(i,j,k)=real(XV(l)*exp(cplx*2*pi*m0*xlon/Lx));
        end
    end

end
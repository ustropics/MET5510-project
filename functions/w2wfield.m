%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: w2wfield.m

% DESCRIPTION: Converts vertical velocity into a 3D wind field over the grid 
% (longitude, latitude, height) for the quasi-geostrophic model.

% INPUT:
% - w: Vertical velocity field
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction

% OUTPUT:
% - wfield: 3D array representing the wind field

% MATH/FUNCTIONS: wfield = w * exp(i * k * x)

% - VARIABLES:
%   - w is the vertical velocity
%   - k is the zonal wavenumber
%   - x is the longitude coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wfield= w2wfield(w,ii,dx)

    global jj kk ll cplx m0 Lx LW
    wfield=zeros(ii+1,jj+1,kk+1);
    
    %% convert linear index to (j,k)
    for l = 1:LW 
        [j,k]=lw2jk(l);
        for i = 1:ii+1
            xlon=(i-1)*dx;
            
            % compute real part of complex exponential for wave structure
            wfield(i,j,k)=real(w(l)*exp(cplx*2*pi*m0*xlon/Lx));
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: XVx2field.m

% Description: Converts the x-derivative of the eigenvector (XVx) into a 
% 3D meridional wind field over the grid in the quasi-geostrophic model, 
% accounting for geostrophic balance.

% INPUT:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction

% OUTPUT:
% - field: 3D array representing the meridional wind field

% MATH/FUNCTIONS: 
% - v' = (g/f₀) ∂ψ/∂x

% VARIABLES:
% - ψ is streamfunction
% - g is gravity
% - f₀ is Coriolis parameter
% - ∂/∂x is derived from XVx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field= XVx2field(XV,ii,dx) 

    global jj kk ll cplx m0 Lx 

    field=zeros(ii+1,jj+1,kk+1); % initialize 3D field

    for l = 1:ll
        [j,k]=l2jk(l); % convert linear index to (j,k)
        for i = 1:ii
            xlon=(i-1)*dx; % longitude coordinate
    
            % Compute the x-derivative using spectral method
            field(i,j,k)=real( (cplx*2*pi*m0/Lx)*XV(l) ...
                *exp(cplx*2*pi*m0*xlon/Lx) );
        
        end
    end
    
end
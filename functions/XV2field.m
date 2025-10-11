%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XV2field.m

% Description: Transforms the 1D eigenvector (XV) into a 3D field 
% (e.g., geopotential height) over longitude, latitude, and height grids using 
% spectral or finite difference methods in the quasi-geostrophic model 

% Math/functions: Field = XV2field(XV, ii, dx) * f₀/g, where 
% XV is the eigenvector
% ii is longitude points
% dx is grid spacing
% f₀ is Coriolis parameter
% and g is gravity

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
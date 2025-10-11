%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: w2wfield.m

% Description: Converts velocity into a 3D wind field over the grid 
% (longitude, latitude, height) for the quasi-geostrophic model

% Math/functions: w_field = w * [u, v, w], where 
% u, v, w are wind components derived from geostrophic balance 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wfield= w2wfield(w,ii,dx)
global jj kk ll cplx m0 Lx LW
wfield=zeros(ii+1,jj+1,kk+1);

% convert linear index to (j,k)
for l = 1:LW 
    [j,k]=lw2jk(l);
    for i = 1:ii+1
        xlon=(i-1)*dx;
        % compute real part of complex exponential for wave structure
        wfield(i,j,k)=real(w(l)*exp(cplx*2*pi*m0*xlon/Lx));
    end
end
end   
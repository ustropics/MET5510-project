%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XVz2field.m

% Description: Converts the z-derivative of the eigenvector (XVz) into a 
% 3D vertical velocity or temperature perturbation field over the grid in 
% the quasi-geostrophic model.

% Math/functions: w' or T' = (f₀ * H / R) ∂ψ/∂z, where 
% ψ is streamfunction
% f₀ is Coriolis parameter
% H is scale heigh
%  R is gas constant
% ∂/∂z is derived from XVz 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field= XVz2field(XV,ii,dx) 
global jj kk ll cplx m0 Lx dz
field=zeros(ii+1,jj+1,kk+1); % initalize 3D field

% convert linear index to (j,k)
for l = 1:ll
    [j,k]=l2jk(l);
    for i = 1:ii
        xlon=(i-1)*dx; % longitude coordinate
        if (k == 1)
            % difference at bottom boundary
            field(i,j,k)=real((XV(jk2l(j,k+1))-XV(l))*exp(cplx*2*pi*m0*xlon/Lx))/dz;
        elseif (k == kk+1)
            % backward difference at top boundary
            field(i,j,k)=real((XV(l)-XV(jk2l(j,k-1)))*exp(cplx*2*pi*m0*xlon/Lx))/dz;
        else
            % inteiro difference 
            field(i,j,k)=real((XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))*exp(cplx*2*pi*m0*xlon/Lx))/(2*dz);
        end
      
    end
end
end  
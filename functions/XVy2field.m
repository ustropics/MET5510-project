%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: XVy2field.m

% Description: Converts the y-derivative of the eigenvector (XVy) into a 
% 3D zonal wind field over the grid in the quasi-geostrophic model, 
% accounting for geostrophic balance.

% Input:
% - XV: Eigenvector representing the streamfunction
% - ii: Number of grid points in the x-direction
% - dx: Grid spacing in the x-direction

% Output:
% - field: 3D array representing the zonal wind field

% Math/functions: u' = -(g/f₀) ∂ψ/∂y

% - Variables:
%   - ψ is streamfunction
%   - g is gravity
%   - f₀ is Coriolis parameter
%   - ∂/∂y is derived from XVy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field= XVy2field(XV,ii,dx) 
    global jj kk ll cplx m0 Lx dy
    
    field=zeros(ii+1,jj+1,kk+1); % initialize 3D field
    
    for l = 1:ll
        [j,k]=l2jk(l);
        for i = 1:ii
            xlon=(i-1)*dx;
            
            %% Apply boundary conditions and compute finite difference
            % For forward difference at southern boundary
            if (j == 2)
                lnh=jk2l(j+1,k);
                field(i,j,k)=-real((XV(lnh)-0.0)*exp(cplx*2*pi*m0*xlon/Lx))/(2*dy); % forward difference at southern boundary

            % For backward difference at northern boundary
            elseif (j == jj)
                lsh=jk2l(j-1,k);
                field(i,j,k)=-real( (0-XV(lsh))*exp(cplx*2*pi*m0*xlon/Lx) )/(2*dy);

            % For centered difference at interior points
            else
                lnh=jk2l(j+1,k);
                lsh=jk2l(j-1,k);
                field(i,j,k)=-real((XV(lnh)-XV(lsh))*exp(cplx*2*pi*m0*xlon/Lx))/(2*dy);
            end
        end
    end
    
    %% Boundary conditions for j=1 and j=jj+1
    for k = 1:kk+1
        lj2=jk2l(2,k); % index at j=2
        ljj=jk2l(jj,k); % index at j=jj
        for i = 1:ii
           xlon=(i-1)*dx;
           field(i,1,k)=-real((XV(lj2)-0.0)*exp(cplx*2*pi*m0*xlon/Lx))/dy;
           field(i,jj+1,k)=-real((0-XV(ljj))*exp(cplx*2*pi*m0*xlon/Lx))/dy;
        end
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: compute_forcing_F1.m

% DESCRIPTION: Vertical shear term forcing – handles all latitude rows
% (boundaries + interior) with the special stencils required at walls.

% INPUT:
%   XV   – eigenvector
%   Ubar – background zonal wind (jj+1 × kk+1)

% OUTPUT:
%   F1 – vector of size LW

function F1 = compute_forcing_F1(XV,Ubar)
    global jj kk cplx m0 Lx dy dz

    LW = (jj+1)*(kk-1);
    F1 = zeros(LW,1);

    % ----- j = 1  (southern boundary)
    j = 1;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1))*...
                (XV(jk2l(3,k))-2*XV(jk2l(2,k)))/(2*dz*dy*dy);
    end

    % ----- j = 2  (first interior point next to southern wall)
    j = 2;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1))*...
                ( -(2*pi*m0/Lx)^2*XV(jk2l(j,k)) + ...
                  (XV(jk2l(3,k))-2*XV(jk2l(2,k)))/dy/dy )/(2*dz);
    end

    % ----- j = jj (first interior point next to northern wall)
    j = jj;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1))*...
                ( -(2*pi*m0/Lx)^2*XV(jk2l(j,k)) + ...
                  (XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/dy/dy )/(2*dz);
    end

    % ----- j = jj+1 (northern boundary)
    j = jj+1;
    for k = 2:kk
        l = jk2lw(j,k);
        F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1))*...
                (XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/(2*dz*dy*dy);
    end

    % ----- interior latitudes 3 ≤ j ≤ jj-1
    for j = 3:jj-1
        for k = 2:kk
            l = jk2lw(j,k);
            F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1))*...
                    ( -(2*pi*m0/Lx)^2*XV(jk2l(j,k)) + ...
                      (XV(jk2l(j+1,k))-2*XV(jk2l(j,k))+XV(jk2l(j-1,k)))/dy/dy )/(2*dz);
        end
    end
end
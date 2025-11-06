 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: compute_forcing_F2.m

% DESCRIPTION: Horizontal shear term forcing (interior only).

% INPUT:  XV – eigenvector
% OUTPUT: F2 – vector of size LW

function F2 = compute_forcing_F2(XV,Ubar)
    global jj kk cplx m0 Lx dy dz

    F2 = zeros((jj+1)*(kk-1),1);

    for j = 2:jj
        for k = 2:kk
            l = jk2lw(j,k);
            F2(l) = -2*(2*pi*m0/Lx)*cplx*...
                    (Ubar(j+1,k)-2*Ubar(j,k)+Ubar(j-1,k))*...
                    (XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz*dy*dy);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: compute_forcing_F3.m

% DESCRIPTION: β-term forcing (interior latitudes only).

% INPUT:  XV – eigenvector (streamfunction)
% OUTPUT: F3 – vector of size LW

function F3 = compute_forcing_F3(XV)
    global jj kk cplx m0 Lx beta dz

    F3 = zeros((jj+1)*(kk-1),1);

    for j = 2:jj
        for k = 2:kk
            l = jk2lw(j,k);
            F3(l) = (2*pi*m0/Lx)*cplx*beta*...
                    (XV(jk2l(j,k+1)) - XV(jk2l(j,k-1)))/(2*dz);
        end
    end
end
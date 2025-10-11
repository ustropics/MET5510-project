% functions/XVz2field.m
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
function wfield= w2wfield(w,ii,dx)
% Convert vertical velocity vector to 3D field
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
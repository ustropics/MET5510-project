% functions/w2ellipse.m
function EW = w2ellipse (w)
global jj kk LW NN2 m0 f0 dy dz Lx 
EW=zeros(LW,1);

% boundary condition at j = 1
j = 1;
for k = 2:kk
    l = jk2lw(j,k);
    ln3=jk2lw(3,k);
    ln2=jk2lw(2,k);
    if (k == 2)
        wdn = 0; % bottom boundary
        wup=w(jk2lw(j,k+1)); % upper level
    elseif (k == kk)
        wdn=w(jk2lw(j,k-1)); % lower level
        wup=0; % top boundary condition
    else 
        wdn=w(jk2lw(j,k-1)); % lower level
        wup=w(jk2lw(j,k+1)); % upper level
    end

    % elliptical operator: includes zonal wavenumber, meridional
    % and vertical derivatives
    EW(l) = -(2*pi*m0/Lx)^2*w(l) ...
            +(w(ln3)-2*w(ln2)+w(l))/dy^2 ...
    +(f0/dz)^2*( wup-2*w(l)+wdn )/NN2;
end

% boundary condition at j = jj + 1
j = jj + 1;
for k = 2:kk
    l = jk2lw(j,k);
    ls3=jk2lw(jj-1,k);
    ls2=jk2lw(jj,k);
    if (k == 2)
        wdn = 0;
        wup=w(jk2lw(j,k+1));
    elseif (k == kk)
        wdn=w(jk2lw(j,k-1));
        wup=0;
    
    else
        wdn=w(jk2lw(j,k-1));
        wup=w(jk2lw(j,k+1));
    end

    EW(l) = -(2*pi*m0/Lx)^2*w(l) ...
            +(w(l)-2*w(ls2)+w(ls3))/dy^2 ...
    +(f0/dz)^2*( wup-2*w(l)+wdn )/NN2;
end

% interior points
for j = 2:jj
    for k = 2:kk
    l = jk2lw(j,k);
    ln1=jk2lw(j+1,k);
    ls1=jk2lw(j-1,k);
    if (k == 2)
        wdn = 0;
        wup=w(jk2lw(j,k+1));
    elseif (k == kk)
        wdn=w(jk2lw(j,k-1));
        wup=0;
    
    else
        wdn=w(jk2lw(j,k-1));
        wup=w(jk2lw(j,k+1));
    end

    EW(l) = -(2*pi*m0/Lx)^2*w(l) ...
            +(w(ln1)-2*w(l)+w(ls1))/dy^2 ...
    +(f0/dz)^2*( wup-2*w(l)+wdn )/NN2;   
    end
end

end
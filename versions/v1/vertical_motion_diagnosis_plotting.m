%from X-vector to streamfunction at grids

load('Rossby_wave_2.mat');
%load('Rossby_wave_6.mat');
tic

global jj kk ll LW BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% grid parameters
ii = 360; % longitude grid points
dx = Lx/ii; % longitude grid spacing

LW=(jj+1)*(kk-1); % matrix for vertical motion

% coordinate arrays
xx=0.0:360/ii:360; % longitude grid
yy=linspace(45-25,45+25,jj+1); % latitude grid
zz=linspace(0.0,10,kk+1); % height grid

% initalize x-vector and compute related field
XV=zeros(ll,1); % initalize x-vector
XV(:)=eigVec3(:,7); % extract 7th eigenvector
QV = B*XV; % PV vector x matrix B

% compute geopotential height field
gpt_h = XV2field(XV,ii,dx)*f0/gg;
[valuemax,indexmax]=max(gpt_h(:));
XV=(10/valuemax)*XV;

%% recalculate fields after normalization
% From Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% Page 7
gpt_h = XV2field(XV,ii,dx)*f0/gg; % geopotential height
pv= XV2field(QV,ii,dx); % PV 
temp = (f0*HH/287)*XVz2field(XV,ii,dx); % temperature

% From Dr. Cai's EigenValue_elementary_analysis_linear_QG_model.pdf
% Pages 9-11 (3D fields for v', u', T')
ug=XVy2field(XV,ii,dx); % zonal wind (u)
vg=XVx2field(XV,ii,dx); % meridional wind (v)

%% initalize G matrix for vertical motion
G=zeros(LW,LW);
for l0 = 1:LW
    w=zeros(LW,1);
    w(l0)=1;
    
    % G matrix
    EW = w2ellipse(w);
    G(:,l0)=EW(:);
end

%% calculating F1, F2, and F3
% from Dr. Cai's Vertical_motion_analysis_linear_QG_model.pdf (slide 1)
F1=zeros(LW,1);
F2=zeros(LW,1);
F3=zeros(LW,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F1, F2, F3 FOR LOOPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute F2 and F3 for interior latitude points (j = 2 to jj)
for j = 2:jj
    for k = 2:kk
    l = jk2lw(j,k);
    F3(l)=(2*pi*m0/Lx)*cplx*beta*(XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz);
    F2(l)=-2*(2*pi*m0/Lx)*cplx*(Ubar(j+1,k)-2*Ubar(j,k)+Ubar(j-1,k)) ...
            *(XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz*dy*dy);
    end
end

%% compute F1 for boundary latitude (j = 1)
j = 1;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(XV(jk2l(3,k))-2*XV(jk2l(2,k)))/(2*dz*dy*dy);
end

%% compute F1 for latitude j = 2
j = 2;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
        (XV(jk2l(3,k))-2*XV(jk2l(2,k)))/dy/dy)/(2*dz);
end

%% compute F1 for latitude j = jj
j = jj ;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
         (XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/dy/dy)/(2*dz);
end

%% compute F1 for latitude j = jj + 1
j = jj + 1;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
          *(XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/(2*dz*dy*dy);
end

%% compute F1 for interior latitudes (j = 3 to jj-1)
for j = 3:jj-1
    for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
         (XV(jk2l(j+1,k))-2*XV(jk2l(j,k))+XV(jk2l(j-1,k)))/dy/dy)/(2*dz);
    end
end

%% vertical velocity (w) using the elliptic equation: G*w = F1 + F2 + F3
w=(f0/NN2)*(G^-1)*(F1+F2+F3);
wfield=w2wfield(w,ii,dx); % Vertical_motion_analysis_linear_QG_model.pdf (page 9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
clear title

% Plot 1: Cross-section of geopotential height at latitude 45N
figure('units','inch','position',[1,1,24,16]);
subplot(2,3,1);
contourf(xx,zz,(squeeze(gpt_h(:,jj/2+1,:)))',-50:2:50,'linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca,'xtick',0:30:360);
set(gca,'ytick',0:1:10);
title('Cross section of geop. height at lat = 45')
set(gca,'Fontsize',16,'Fontweight','Bold');
% 

% Plot 2: Cross-section of vertical motion at latitude 45N
subplot(2,3,4);
ttt=squeeze(wfield(:,jj/2+1,:));
ttt(abs(ttt)<1.0e-8)=0; %set all noise data to zero
contourf(xx,zz,ttt','linestyle', 'none');
%contourf(xx,zz,ttt',-50:0.2:50,'linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca,'xtick',0:30:360);
set(gca,'ytick',0:1:10);
title('Cross section of vertical motions at lat = 45')
set(gca,'Fontsize',16,'Fontweight','Bold');

% Plot 3: Temperature and geopotential height at 300 hPa (~7 km)
subplot(2,3,2)
ttt=squeeze(temp(:,:,kk/2+1+10));
ttt(abs(ttt)<1.0e-8)=0; %set all noise data to zero
contourf(xx,yy,ttt',-50:0.2:50,'linestyle', 'none');
hold on;
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1+10)))',0:2:50,'w');
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1+10)))',-50:2:-2,'--w');
colorbar;
set(gca,'clim',[-1 1]);
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Temp/height (shading/contour) at 300 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

% Plot 4: Temperature and geopotential height at 700 hPa (~3 km)
subplot(2,3,5)
ttt=squeeze(temp(:,:,kk/2+1-10));
ttt(abs(ttt)<1.0e-8)=0; %set all noise data to zero
contourf(xx,yy,ttt',-50:0.2:50,'linestyle', 'none');
hold on;
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1-10)))',0:2:50,'w');
contour(xx,yy,(squeeze(gpt_h(:,:,kk/2+1-10)))',-50:2:-2,'--w');
colorbar;
set(gca,'clim',[-1 1]);
%set(gca,'clim',[-max(ttt(:))-0.000001 max(ttt(:))+0.000001]);
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Temp/height (shading/contour) at 700 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

% Plot 5: Vertical motion at 300 hPa
subplot(2,3,3)
contourf(xx,yy,(squeeze(wfield(:,:,kk/2+1+10)))','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vertical motion at 300 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

% Plot 6: Vertical motion at 800 hPa (~2 km)
subplot(2,3,6)
contourf(xx,yy,(squeeze(wfield(:,:,kk/2+1-10)))','linestyle', 'none');
% contourf(xx,yy,(squeeze(vg(:,:,kk/2+1)))',-10:0.1:10);
colorbar;
xlabel('Longitude')
ylabel('latitude')
%set(gca,'xlim',[40 110])
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vetical motion at 800 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function jk2lw: Function to convert (j,k) indices to linear index for vertical motion
function l = jk2lw(j,k)
global jj  

l = j + (k-2)*(jj+1);

end

%% Function lw2jk: Function to convert linear index to (j,k) 
% indices for vertical motion for k floor and j (jj + 1)
function [j,k] = lw2jk(l)
global jj

k = floor((l-1)/(jj+1))+2;
j = l - (k-2)*(jj+1);
end

%% Function jk2l: to convert (j,k) indices to linear index for streamfunction
% for l
function l = jk2l(j,k)
global jj  

l = j-1 + (k-1)*(jj-1);

end

%% Function to convert linear index to (j,k) indices for streamfunction
% for vertical motion for k floor and j (jj - 1)
function [j,k] = l2jk(l)
global jj
k = floor((l-1)/(jj-1))+1;
j = l+1 - (k-1)*(jj-1);
end

%% Function: EW ellipse array for vertical motion
% From Dr. Cai's Vertical_motion_analysis_linear_QG_model.pdf
% slides 3-4
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

%% Function wfield: convert vertical velocity vector to 3D field
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

%$ Function XV2field: covert X-vector to 3D field (streamfunction/GPH)
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

%% Function Xvz2field: compute vertical derivative of x-vector field
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

%% Function XVx2field: compute meridional wind (v) from X-vector
function field= XVx2field(XV,ii,dx) 

global jj kk ll cplx m0 Lx 

field=zeros(ii+1,jj+1,kk+1); % initialize 3D field

for l = 1:ll
    [j,k]=l2jk(l); % convert linear index to (j,k)
    for i = 1:ii
        xlon=(i-1)*dx; % longitude coordinate
 
        field(i,j,k)=real( (cplx*2*pi*m0/Lx)*XV(l) ...
            *exp(cplx*2*pi*m0*xlon/Lx) );
      
    end
end
end  

%% Function XVy2field: compute zonal wind (u) from X-vector
function field= XVy2field(XV,ii,dx) 
global jj kk ll cplx m0 Lx dy

field=zeros(ii+1,jj+1,kk+1); % initialize 3D field

for l = 1:ll
    [j,k]=l2jk(l);
    for i = 1:ii
        xlon=(i-1)*dx;
        
        if (j == 2)
            lnh=jk2l(j+1,k);
            field(i,j,k)=-real((XV(lnh)-0.0)*exp(cplx*2*pi*m0*xlon/Lx))/(2*dy);
        elseif (j == jj)
            lsh=jk2l(j-1,k);
            field(i,j,k)=-real( (0-XV(lsh))*exp(cplx*2*pi*m0*xlon/Lx) )/(2*dy);
        else
            lnh=jk2l(j+1,k);
            lsh=jk2l(j-1,k);
            field(i,j,k)=-real((XV(lnh)-XV(lsh))*exp(cplx*2*pi*m0*xlon/Lx))/(2*dy);
        end
        
    end
end

% Boundary conditions for j=1 and j=jj+1
for k = 1:kk+1
    lj2=jk2l(2,k);
    ljj=jk2l(jj,k);
    for i = 1:ii
       xlon=(i-1)*dx;
       field(i,1,k)=-real((XV(lj2)-0.0)*exp(cplx*2*pi*m0*xlon/Lx))/dy;
       field(i,jj+1,k)=-real((0-XV(ljj))*exp(cplx*2*pi*m0*xlon/Lx))/dy;
    end

end
end  
    





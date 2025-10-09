% main_script.m
% Load data and set up global variables
load('data/Rossby_wave_2.mat');
%load('Rossby_wave_6.mat');
tic

global jj kk ll LW BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx HH gg

% Add functions directory to MATLAB path
addpath('functions');

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

% initialize x-vector and compute related field
XV=zeros(ll,1); % initialize x-vector
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

%% initialize G matrix for vertical motion
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

% Plot 2: Cross-section of vertical motion at latitude 45N
subplot(2,3,4);
ttt=squeeze(wfield(:,jj/2+1,:));
ttt(abs(ttt)<1.0e-8)=0; %set all noise data to zero
contourf(xx,zz,ttt','linestyle', 'none');
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
colorbar;
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vertical motion at 800 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');
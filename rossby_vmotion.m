%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_vmotion.m

% Description: Script for loading Eady wave data from 'eady_wave.mat' and 
% computing vertical motion diagnostics in the quasi-geostrophic model, 
% including geopotential height, temperature, wind components, and vertical 
% velocity fields, with subsequent plotting of cross-sections and pressure 
% level diagnostics.

% Functions used: 
% - jk2lw: Converts 2D indices (j, k) to a weighted linear index.
% - XV2field: Transforms eigenvector to 3D field.
% - XVz2field: Converts z-derivative of eigenvector to 3D field.
% - XVy2field: Converts y-derivative of eigenvector to 3D field.
% - XVx2field: Converts x-derivative of eigenvector to 3D field.
% - w2ellipse: Transforms weight vector to an elliptical representation.
% - w2wfield: Converts weight vector to a 3D wind field.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Rossby_wave_6.mat');
tic

global jj kk ll LW BPVy NN2 f0 dy dz m0 Lx Ubar f0 beta cplx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Grid parameters
ii = 360; % longitude grid points
dx = Lx/ii; % longitude grid spacing

LW=(jj+1)*(kk-1); % matrix for vertical motion

%% Coordinate arrays
xx=0.0:360/ii:360; % longitude grid
yy=linspace(45-25,45+25,jj+1); % latitude grid
zz=linspace(0.0,10,kk+1); % height grid

%% Initalize x-vector and compute related field
XV=zeros(ll,1); % initalize x-vector
XV(:)=eigVec3(:,7); % extract 7th eigenvector
QV = B*XV; % PV vector x matrix B

% compute geopotential height field
gpt_h = XV2field(XV,ii,dx)*f0/gg;
[valuemax,indexmax]=max(gpt_h(:));
XV=(10/valuemax)*XV;

%% Recalculate fields after normalization
% From EigenValue_elementary_analysis_linear_QG_model.pdf
% Page 7
gpt_h = XV2field(XV,ii,dx)*f0/gg; % geopotential height
pv= XV2field(QV,ii,dx); % PV 
temp = (f0*HH/287)*XVz2field(XV,ii,dx); % temperature

% From EigenValue_elementary_analysis_linear_QG_model.pdf
% Pages 9-11 (3D fields for v', u', T')
ug=XVy2field(XV,ii,dx); % zonal wind (u)
vg=XVx2field(XV,ii,dx); % meridional wind (v)

%% Initalize G matrix for vertical motion
G=zeros(LW,LW);
for l0 = 1:LW
    w=zeros(LW,1);
    w(l0)=1;
    
    % G matrix
    EW = w2ellipse(w);
    G(:,l0)=EW(:);
end

%% Calculating F1, F2, and F3
% from Vertical_motion_analysis_linear_QG_model.pdf (slide 1)
F1=zeros(LW,1);
F2=zeros(LW,1);
F3=zeros(LW,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% F1, F2, F3 FOR LOOPS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute F2 and F3 for interior latitude points (j = 2 to jj)
for j = 2:jj
    for k = 2:kk
    l = jk2lw(j,k);
    F3(l)=(2*pi*m0/Lx)*cplx*beta*(XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz);
    F2(l)=-2*(2*pi*m0/Lx)*cplx*(Ubar(j+1,k)-2*Ubar(j,k)+Ubar(j-1,k)) ...
            *(XV(jk2l(j,k+1))-XV(jk2l(j,k-1)))/(2*dz*dy*dy);
    end
end

%% Compute F1 for boundary latitude (j = 1)
j = 1;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(XV(jk2l(3,k))-2*XV(jk2l(2,k)))/(2*dz*dy*dy);
end

%% Compute F1 for latitude j = 2
j = 2;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
        (XV(jk2l(3,k))-2*XV(jk2l(2,k)))/dy/dy)/(2*dz);
end

%% Compute F1 for latitude j = jj
j = jj ;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
         (XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/dy/dy)/(2*dz);
end

%% Compute F1 for latitude j = jj + 1
j = jj + 1;
for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
          *(XV(jk2l(jj-1,k))-2*XV(jk2l(jj,k)))/(2*dz*dy*dy);
end

%% Compute F1 for interior latitudes (j = 3 to jj-1)
for j = 3:jj-1
    for k = 2:kk
    l = jk2lw(j,k);
    F1(l) = 2*(2*pi*m0/Lx)*cplx*(Ubar(j,k+1)-Ubar(j,k-1)) ...
        *(-(2*pi*m0/Lx)^2*XV(jk2l(j,k))+ ...
         (XV(jk2l(j+1,k))-2*XV(jk2l(j,k))+XV(jk2l(j-1,k)))/dy/dy)/(2*dz);
    end
end

%% Vertical velocity (w) using the elliptic equation: G*w = F1 + F2 + F3
w=(f0/NN2)*(G^-1)*(F1+F2+F3);
wfield=w2wfield(w,ii,dx); % Vertical_motion_analysis_linear_QG_model.pdf (page 9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
clear title % fix annoying bug from rossby_plot.m

% Create plots folder if it doesn't exist
if ~exist('plots', 'dir')
    mkdir('plots');
end

%% Plot 1-6: Combined figure with subplots
fig = figure('units','inch','position',[1,1,24,16], 'Visible', 'off');

subplot(2,3,1);
contourf(xx,zz,(squeeze(gpt_h(:,jj/2+1,:)))',-50:2:50,'linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca,'xtick',0:30:360);
set(gca,'ytick',0:1:10);
title('Cross section of geop. height at lat = 45')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,4);
ttt=squeeze(wfield(:,jj/2+1,:));
ttt(abs(ttt)<1.0e-8)=0; % set all noise data to zero
contourf(xx,zz,ttt','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Height')
set(gca,'xtick',0:30:360);
set(gca,'ytick',0:1:10);
title('Cross section of vertical motions at lat = 45')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,2)
ttt=squeeze(temp(:,:,kk/2+1+10));
ttt(abs(ttt)<1.0e-8)=0;
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

subplot(2,3,5)
ttt=squeeze(temp(:,:,kk/2+1-10));
ttt(abs(ttt)<1.0e-8)=0; 
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

subplot(2,3,3)
contourf(xx,yy,(squeeze(wfield(:,:,kk/2+1+10)))','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vertical motion at 300 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

subplot(2,3,6)
contourf(xx,yy,(squeeze(wfield(:,:,kk/2+1-10)))','linestyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('latitude')
set(gca,'xtick',0:60:360);
set(gca,'ytick',25:5:65);
title('Vetical motion at 800 hPa')
set(gca,'Fontsize',16,'Fontweight','Bold');

% Save the figure
saveas(fig, 'plots/vertical_motion_diagnosis.png');
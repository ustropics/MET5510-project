%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILENAME: calc.m
% DESCRIPTION: Post-processing driver for HWME model – now processes 
%              modes 1 to 12 automatically (no need to set n_mode in cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; close all; clc;
global jj kk ll xx yy zz HH ii BPVy NN2 f0 gg dx dy m0 dz Lx Theta0 ...
    Ubar beta cplx hlat hlevel time LW fig_path

tic

addpath(['functions', filesep])

%% Load parameters (cfg no longer needs n_mode)
params = cfg();
jj = params.jj; kk = params.kk; ll = params.ll; Lx = params.Lx; LW = params.LW;
f0 = params.f0; beta = params.beta; gg = params.gg; Theta0 = params.Theta0;
HH = params.HH; dy = params.dy; dz = params.dz; NN2 = params.NN2;
cplx = params.cplx; m0 = params.m0; fig_path = params.plot_dir;
ii = params.ii; dx = params.dx; 
xx = params.xx; yy = params.yy; zz = params.zz;
hlat = params.hlat; time = params.time;

%% Load eigenmode data (once)
data_file = fullfile(params.data_dir, params.data_filename);
load(data_file, 'eigVec2', 'eigVal2', 'B', 'Ubar', 'BPVy');

fprintf('Loaded eigenmodes from: %s\n', data_file);
fprintf('Starting calculation for modes 1 to 12...\n\n');

%% Build the inversion matrix G ONLY ONCE (this was the source of the error!)
G = Gmatrix();   % Size depends only on vertical grid → mode-independent

%% Ensure output directory exists
if ~exist(params.plot_dir, 'dir')
    mkdir(params.plot_dir);
end

%% =======================================================================
%% MAIN LOOP OVER MODES 1 to 12
%% =======================================================================
for n_mode = 1:12
    
    fprintf('=== Processing mode n = %2d ===\n', n_mode);

    %% Select eigenmode
    XV = eigVec2(:, n_mode);
    omega       = imag(eigVal2(n_mode));
    growth_rate = real(eigVal2(n_mode));
    
    if growth_rate <= 0
        fprintf('  Mode %d is neutral or decaying (σ = %.2e) – skipping.\n', ...
                n_mode, growth_rate);
        continue;
    end

    phase_speed = -omega / (2*pi*m0/Lx);                    % m/s
    eFolding    = (1/growth_rate)/86400;                    % days

    %% Normalise so that maximum |ψ| corresponds to 10 m geopotential height anomaly
    psi_3D = XV2field(XV, ii, dx);                          % streamfunction field
    valuermax = max(abs(psi_3D * f0 / gg), [], 'all');
    
    if valuermax == 0 || isinf(valuermax) || isnan(valuermax)
        fprintf('  Warning: Mode %d has zero amplitude – skipping.\n', n_mode);
        continue;
    end
    
    scaling = 10 / valuermax;
    XV = scaling * XV;                                      % normalised eigenvector

    %% Compute basic fields
    QV  = B * XV;
    XVy = XV2XVy(XV);
    XVx = XV2XVx(XV);
    XVz = XV2XVz(XV);

    gpt_h   =  XV2field(XV , ii, dx) * f0/gg;              % geopot. height anomaly (m)
    temp    = (f0*HH/287) * XV2field(XVz, ii, dx);          % temperature anomaly (K)
    ug      = -XV2field(XVy, ii, dx);                       % geostrophic zonal wind
    vg      =  XV2field(XVx, ii, dx);                       % geostrophic meridional wind
    pvfield =  XV2field(QV , ii, dx);                       % PV anomaly

    %% Omega equation forcing terms (now compatible with pre-computed G)
    F3 = F3matrix(XV);
    F2 = F2matrix(XV, Ubar);
    F1 = F1matrix(XV, Ubar);

    w_vec   = w2vec(G, F1, F2, F3);                         % interior solution
    wfield  = w2field(w_vec, ii, dx);                       % full 3D vertical velocity

    %% Hovmöller diagrams at selected latitude and two levels
    gpt_h_hovmoler1  = XV2streamxtime(XV, ii, dx, omega, hlat, 1)  * f0/gg;
    gpt_h_hovmoler51 = XV2streamxtime(XV, ii, dx, omega, hlat, 51) * f0/gg;
    ug_hovmoler1     = XV2ugxtime(XVy, ii, dx, omega, hlat, 1);
    ug_hovmoler51    = XV2ugxtime(XVy, ii, dx, omega, hlat, 51);

    %% Optional: store un-normalised eigenvector too (some users like both)
    XV_raw = eigVec2(:, n_mode);

    %% Save results – unique filename per mode
    calc_filename = sprintf('calc_wave-m%d_eMode-%02d.mat', m0, n_mode);
    calcFile = fullfile(params.data_dir, calc_filename);

    save(calcFile, ...
        'gpt_h', 'temp', 'ug', 'vg', 'pvfield', 'wfield', ...
        'gpt_h_hovmoler1', 'gpt_h_hovmoler51', ...
        'ug_hovmoler1', 'ug_hovmoler51', ...
        'xx', 'yy', 'zz', 'time', 'Ubar', 'BPVy', ...
        'XV', 'XV_raw', 'm0', 'n_mode', ...
        'growth_rate', 'omega', 'phase_speed', 'eFolding', ...
        'hlat', 'params', '-v7.3');

    fprintf('  Saved: %s\n', calc_filename);
    
end % end of for loop

fprintf('\nAll 12 modes processed successfully!\n');
toc
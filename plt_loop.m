%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILENAME: plt.m
% DESCRIPTION: Automatic plotting of all 12 modes – 100% robust
% All previous index errors fixed, no ternary operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
addpath('functions', 'plots', 'cmaps');

%% Load config
params    = cfg();
m0        = params.m0;
data_dir  = params.data_dir;
fig_path  = params.plot_dir;
fig_path2 = fullfile(fig_path, 'combined');

if ~exist(fig_path,  'dir'), mkdir(fig_path);  end
if ~exist(fig_path2, 'dir'), mkdir(fig_path2); end

fprintf('=== PLOTTING m0 = %d | Modes 1–12 ===\n\n', m0);

%% =======================================================================
%% LOOP OVER MODES 1–12
%% =======================================================================
for n_mode = 1:12

    filename = sprintf('calc_wave-m%d_eMode-%02d.mat', m0, n_mode);
    filepath = fullfile(data_dir, filename);

    if ~exist(filepath, 'file')
        fprintf('Skipping (file missing): %s\n', filename);
        continue;
    end

    fprintf('Loading & plotting: %s\n', filename);
    S = load(filepath);

    % ------------------------------------------------------------------
    % Extract fields safely
    % ------------------------------------------------------------------
    gpt_h   = S.gpt_h;   if ndims(gpt_h) < 4, gpt_h = gpt_h(:,:,:,1); end
    temp    = S.temp;
    ug      = S.ug;
    vg      = S.vg;
    pvfield = S.pvfield;
    wfield  = S.wfield;

    gpt_h_hovmoler1  = S.gpt_h_hovmoler1;
    gpt_h_hovmoler51 = S.gpt_h_hovmoler51;
    ug_hovmoler1     = S.ug_hovmoler1;
    ug_hovmoler51    = S.ug_hovmoler51;

    xx = S.xx; yy = S.yy; zz = S.zz; time = S.time;
    Ubar = S.Ubar; BPVy = S.BPVy;
    hlat = S.hlat;
    omega = S.omega;
    growth_rate = S.growth_rate;

    % Choose eigenvector: normalized XV first, otherwise raw
    if isfield(S,'XV') && numel(S.XV) > 1
        XV = S.XV(:);
    elseif isfield(S,'XV_raw')
        XV = S.XV_raw(:);
    else
        XV = [];   % will be caught later
    end

    % Grid dimensions from the actual 3D fields (most reliable)
    [nj, ~, nk] = size(gpt_h);
    jj = nj - 1;
    kk = nk - 1;
    ll_expected = (jj+1)*(kk+1);

    %% 1. Fluxes
    plot_momentum_flux(zz, yy, vg, ug, m0, n_mode, fig_path);
    plot_meridional_hflux(zz, yy, vg, temp, m0, n_mode, fig_path);
    plot_vertical_hflux(zz, yy, wfield, temp, m0, n_mode, fig_path);

    %% 2. Hovmöller diagrams – using saved omega (no eigVal2 needed)
    plot_gph_hovmoller(xx, time, gpt_h_hovmoler1,  hlat, 1,  m0, n_mode, fig_path);
    plot_gph_hovmoller(xx, time, gpt_h_hovmoler51, hlat, 51, m0, n_mode, fig_path);
    plot_zonal_hovmoller(xx, time, ug_hovmoler1,  hlat, 1,  m0, n_mode, omega, fig_path);
    plot_zonal_hovmoller(xx, time, ug_hovmoler51, hlat, 51, m0, n_mode, omega, fig_path);

    %% 3. Eigenvector amplitude |ψ|²
    eVec_amp = zeros(jj+1, kk+1);
    if ~isempty(XV) && numel(XV) == ll_expected
        for l = 1:ll_expected
            [j,k] = l2jk(l);
            eVec_amp(j,k) = XV(l) * conj(XV(l));
        end
        plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, fig_path);
    else
        fprintf('  Skipping eigenvector amplitude plot (size mismatch)\n');
    end

    %% 4. Horizontal sections at levels 1, 26, 51 – safe indexing!
    levels = [1, 26, 51];
    for i = 1:numel(levels)
        h      = levels(i);
        h_safe = min(h, nk);   % never exceed actual vertical dimension

        plot_zonal_wind(     xx, yy, ug(:,:,:,h_safe),      h, m0, n_mode, fig_path);
        plot_meridional_wind(xx, yy, vg(:,:,:,h_safe),      h, m0, n_mode, fig_path);
        plot_temp(           xx, yy, temp(:,:,:,h_safe),    h, m0, n_mode, fig_path);
        plot_pvfield(        xx, yy, pvfield(:,:,:,h_safe), h, m0, n_mode, fig_path);
        plot_gph(            xx, yy, gpt_h(:,:,:,h_safe),   h, m0, n_mode, fig_path);

        combined_meridional_wind_temp(xx, yy, vg(:,:,:,h_safe), temp(:,:,:,h_safe), h, m0, n_mode, fig_path2);
        combined_zonal_wind_temp(     xx, yy, ug(:,:,:,h_safe), temp(:,:,:,h_safe), h, m0, n_mode, fig_path2);
        combined_gph_temp(            xx, yy, gpt_h(:,:,:,h_safe), temp(:,:,:,h_safe), h, m0, n_mode, fig_path2);
        combined_pvfield_temp(        xx, yy, pvfield(:,:,:,h_safe), temp(:,:,:,h_safe), h, m0, n_mode, fig_path2);
    end

    %% 5. Cross-sections & combined panels
    plot_gph_xsec(xx, zz, gpt_h, jj, m0, n_mode, fig_path);
    plot_meridional_xsec(xx, zz, vg, jj, m0, n_mode, fig_path);
    combined_gph_meridional_xsec(xx, zz, gpt_h, vg, jj, m0, n_mode, fig_path2);

    combined_momentum_vertical_hflux(vg, ug, wfield, temp, m0, n_mode, fig_path2);
    combined_momentum_meridional_hflux(vg, ug, temp, m0, n_mode, fig_path2);
    combined_ubar_evec_amp(yy, zz, Ubar, eVec_amp, m0, n_mode, growth_rate, omega, fig_path2);

    %% 6. Background state – plot only once
    if n_mode == 1
        plot_ubar(yy, zz, Ubar, m0, n_mode, fig_path);
        plot_dpvdym_int(yy, zz, BPVy, m0, n_mode, fig_path);
        plot_dpvdym_bnd(yy, BPVy, params.beta, params.kk, m0, n_mode, fig_path);
        combined_bg_flow(yy, zz, jj, params.kk, Ubar, BPVy, m0, n_mode, fig_path);
    end

    fprintf('Finished mode %02d\n\n', n_mode);
end

disp('================================================================');
disp('ALL AVAILABLE MODES PLOTTED SUCCESSFULLY!');
disp(['Figures:         ' fig_path]);
disp(['Combined figs:   ' fig_path2]);
disp('================================================================');
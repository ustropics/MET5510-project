%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_one
% 
% DESCRIPTION: Run a single plotting function from plt.m
%
% USAGE EXAMPLES:
%   plot_one('meridional_xsec')                 % no level
%   plot_one('zonal_wind', 25)                  % with level
%   plot_one('combined_zonal_wind_temp', 51)
%
% FUNCTION: This loads cfg(), loads the precomputed .mat file, and calls 
% the chosen plotting function with correct arguments.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plt_one(plotname, h)

addpath(['cmaps', filesep])
    % ---------------------------------------------------------------------
    % 1. Add required paths and load config/initial data
    % ---------------------------------------------------------------------

    % If called with no arguments or 'list', show available plots
    if nargin == 0 || (ischar(plotname) && strcmpi(plotname, 'list'))
        listAvailablePlots();
        return;
    end

    addpath(['functions', filesep])
    addpath(['plots', filesep])

    % Load config and precomputed data
    params = cfg();
    load(fullfile(params.data_dir, params.calc_filename));

    % Unpack shared parameters
    xx = params.xx; yy = params.yy; zz = params.zz;
    jj = params.jj; m0 = params.m0; n_mode = params.n_mode;
    time = params.time; hlat = params.hlat;
    fig_path = params.plot_dir;
    fig_path2 = fullfile(fig_path, 'combined');

    % Ensure figure directories exist
    if ~exist(fig_path, 'dir'), mkdir(fig_path); end
    if ~exist(fig_path2, 'dir'), mkdir(fig_path2); end

    % Checks if a height level was given before proceeding
    if nargin < 2
        h = [];
    end

    % ---------------------------------------------------------------------
    % 2. Checks if height level was given, finds the correct plot function
    % ---------------------------------------------------------------------
    switch plotname

        % === Core diagnostic plots ===
        case 'gph_xsec'
            plot_gph_xsec(xx, zz, gpt_h, jj, m0, n_mode, fig_path);

        case 'momentum_flux'
            plot_momentum_flux(zz, yy, vg, ug, m0, n_mode, fig_path);

        case 'meridional_hflux'
            plot_meridional_hflux(zz, yy, vg, temp, m0, n_mode, fig_path);

        case 'vertical_hflux'
            plot_vertical_hflux(zz, yy, wfield, temp, m0, n_mode, fig_path);

        % === Eigenvector amplitude ===
        case 'evec_amp'
            eVec_amp = zeros(jj + 1, params.kk + 1);
            for l = 1 : params.ll
                [j, k] = l2jk(l);
                eVec_amp(j, k) = XV(l) .* conj(XV(l));
            end
            plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, fig_path);

        % === Combined diagnostic plots ===
        case 'combined_momentum_and_vertical_hflux'
            combined_momentum_and_vertical_hflux(vg, ug, wfield, temp, m0, n_mode, fig_path2);

        case 'combined_momentum_and_meridional_hflux'
            combined_momentum_and_meridional_hflux(vg, ug, temp, m0, n_mode, fig_path2);

        case 'combined_gph_meridional_xsec'
            combined_gph_meridional_xsec(xx, zz, gpt_h, vg, jj, m0, n_mode, fig_path2)

        case 'combined_ubar_with_evec_amp'
            eVec_amp = zeros(jj + 1, params.kk + 1);
            for l = 1 : params.ll
                [j, k] = l2jk(l);
                eVec_amp(j, k) = XV(l) .* conj(XV(l));
            end
            combined_ubar_with_evec_amp(yy, zz, Ubar, eVec_amp, m0, n_mode, ...
                                        growth_rate, omega, fig_path2);

        % === Level-dependent plots ===
        case 'zonal_wind'
            requireLevel(h);
            plot_zonal_wind(xx, yy, ug, h, m0, n_mode, fig_path);

        case 'meridional_wind'
            requireLevel(h);
            plot_meridional_wind(xx, yy, vg, h, m0, n_mode, fig_path);

        case 'temperature'
            requireLevel(h);
            plot_temp(xx, yy, temp, h, m0, n_mode, fig_path);

        case 'pvfield'
            requireLevel(h);
            plot_pvfield(xx, yy, pvfield, h, m0, n_mode, fig_path);

        case 'combined_meridional_wind_temp'
            requireLevel(h);
            combined_meridional_wind_temp(xx, yy, vg, temp, h, m0, n_mode, fig_path2);

        case 'combined_zonal_wind_temp'
            requireLevel(h);
            combined_zonal_wind_temp(xx, yy, ug, temp, h, m0, n_mode, fig_path2);

        case 'combined_gph_temp'
            requireLevel(h);
            combined_gph_temp(xx, yy, gpt_h, temp, h, m0, n_mode, fig_path2);
        
        case 'combined_pvfield_temp'
            requireLevel(h);
            combined_pvfield_temp(xx, yy, pvfield, temp, h, m0, n_mode, fig_path2)

        case 'gph_hovmoller'
            requireLevel(h);
            var_name = sprintf('gpt_h_hovmoler%d', h);
            data = eval(var_name);
            plot_gph_hovmoller(xx, time, data, hlat, h, m0, n_mode, fig_path);

        case 'zonal_hovmoller'
            requireLevel(h);
            var_name = sprintf('ug_hovmoler%d', h);
            data = eval(var_name);
            plot_zonal_hovmoller(xx, time, data, hlat, h, m0, n_mode, fig_path);    

        % === Top boundary and background ===
        case 'gph'
            requireLevel(h);
            plot_gph(xx, yy, gpt_h, h, m0, n_mode, fig_path);

        case 'meridional_xsec'
            plot_meridional_xsec(xx, zz, vg, jj, m0, n_mode, fig_path);

        case 'ubar'
            plot_ubar(yy, zz, Ubar, m0, n_mode, fig_path);

        case 'dpvdym_int'
            plot_dpvdym_int(yy, zz, BPVy, m0, n_mode, fig_path);

        case 'dpvdym_bnd'
            plot_dpvdym_bnd(yy, BPVy, params.beta, params.kk, m0, n_mode, fig_path);

        case 'combined_bg_flow'
            combined_bg_flow(yy, zz, jj, params.kk, Ubar, BPVy, m0, n_mode, fig_path);

        % === Plot was not found ===    
        otherwise
            fprintf('\n Unknown plot name: "%s", here is a list of plots \n', plotname);
            listAvailablePlots();
            return;
    end

    % === Plot succeeded ===
    fprintf('\nPlot "%s" completed successfully.\n', plotname);
end

% -------------------------------------------------------------------------
% 3. Helper functions
% -------------------------------------------------------------------------

% Plot requires h argument
function requireLevel(h)
    if isempty(h)
        error('This plot requires a level argument h, e.g. plot_one(''zonal_wind'', 25)');
    end
end

% List all available plot names
function names = listAvailablePlots()
    names = {
        'gph_xsec', 'momentum_flux', 'meridional_hflux', 'vertical_hflux', ...
        'gph_hovmoller', 'zonal_hovmoller', ...
        'evec_amp', 'combined_momentum_vertical_hflux', ...
        'combined_momentum_meridional_hflux', 'combined_ubar_evec_amp', ...
        'zonal_wind', 'meridional_wind', 'temperature', ...
        'combined_meridional_wind_temp', 'combined_zonal_wind_temp', ...
        'gph', 'pvfield', 'meridional_xsec', 'ubar', ...
        'dpvdym_int', 'dpvdym_bnd', 'combined_bg_flow', ...
        'combined_gph_temp','combined_pvfield_temp', 'combined_gph_meridional_xsec'
    };

    names = sort(names); % sort them so they're alphabetical
    fprintf('Available plot names:\n');
    
    % loop over names
    for i = 1:numel(names)
        fprintf('  %2d. %s\n', i, names{i});
    end

    fprintf('\nUse plot_one(''name'') or plot_one(''name'', level)\n\n');
end


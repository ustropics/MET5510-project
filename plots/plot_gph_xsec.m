%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_gph_xsec.m

% DESCRIPTION: Plots geopotential height perturbation at mid-latitude
%              (jj/2+1) across longitude and height using filled contours.

% INPUT:
%   xx       - Longitude coordinates (degrees), size (ii+1)
%   zz       - Height coordinates (km), size (kk+1)
%   gpt_h    - Geopotential height field (ii+1 x jj+1 x kk+1 array, m)
%   jj       - Number of latitude grid points
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

% OUTPUT:
%   Saves: fig_path/gph_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_gph_xsec(xx, zz, gpt_h, jj, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract mid-latitude slice and compute limits
    %% --------------------------------------------------------------------
    
    % Get 2d data for geopotential height contour at mid-latitude
    data = squeeze(gpt_h(:, floor(jj/2) + 1, :));  % (ii+1 x kk+1)
    vmin = min(data(:));
    vmax = max(data(:));
    
    fprintf('\nGeopotential height contour at mid-latitude:\n')
    fprintf('Geopotential height (mid-lat) - Max: %.0f, Min: %.of [m]\n', vmax, vmin);
    
    step = 0.2;
    vmin = floor(vmin / step) * step;
    vmax = ceil(vmax / step) * step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4,2,16,12],'Visible','off');
    contourf(xx, zz, data', 'LineStyle','none');
    hold on
    contour(xx, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off
    
    colormap(cmap_coolwarm(256));
    colorbar;
    caxis([vmin vmax]);

    xlabel('Longitude (degrees)');
    ylabel('Height (km)');
    set(gca, 'XTick', 0:30:360, 'YTick', 0:2:10);

    title_str = ['Geopotential Height Perturbation (Mid-Latitude)' newline ...
                 'zonal wave # = ', num2str(m0), ', eMode # = ', num2str(n_mode)];
    title(title_str);

    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ['gph_xsec', ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving geopotential height plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end
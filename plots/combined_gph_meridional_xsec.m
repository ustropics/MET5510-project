%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_gph_meridional_xsec.m
%
% DESCRIPTION: Plots geopotential height perturbation (filled contours) at
%              mid-latitude with meridional wind (black contour lines)
%              overlaid.
%
% INPUT:
%   xx       - Longitude coordinates (degrees), size (ii+1)
%   zz       - Height coordinates (km), size (kk+1)
%   gpt_h    - Geopotential height field (ii+1 x jj+1 x kk+1 array, m)
%   vg       - Meridional wind field (ii+1 x jj+1 x kk+1 array, m/s)
%   jj       - Number of latitude grid points
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure
%
% OUTPUT:
%   Saves: fig_path/combined_gph_meridional_xsec_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_gph_meridional_xsec(xx, zz, gpt_h, vg, jj, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract mid-latitude slices
    %% --------------------------------------------------------------------
    lat_idx = floor(jj/2) + 1;

    gph_data = squeeze(gpt_h(:, lat_idx, :));  % (ii+1 x kk+1)
    vg_data  = squeeze(vg(:, lat_idx, :));     % (ii+1 x kk+1)

    %% --------------------------------------------------------------------
    %% 2. Compute color limits for GPH (filled)
    %% --------------------------------------------------------------------
    vmin_gph = min(gph_data(:));
    vmax_gph = max(gph_data(:));
    
    step_gph = 0.2;
    vmin_gph = floor(vmin_gph / step_gph) * step_gph;
    vmax_gph = ceil(vmax_gph / step_gph) * step_gph;

    fprintf('\nGeopotential Height (mid-lat): Max = %.1f, Min = %.1f [m]\n', vmax_gph, vmin_gph);

    %% --------------------------------------------------------------------
    %% 3. Compute contour levels for meridional wind (lines)
    %% --------------------------------------------------------------------
    vmin_vg = min(vg_data(:));
    vmax_vg = max(vg_data(:));
    
    step_vg = 0.05;
    vmin_vg = floor(vmin_vg / step_vg) * step_vg;
    vmax_vg = ceil(vmax_vg / step_vg) * step_vg;

    % Generate contour levels (avoid zero if possible unless meaningful)
    vg_levels = vmin_vg:step_vg*2:vmax_vg;  % wider spacing for lines
    if ~any(vg_levels == 0) && (vmin_vg < 0 && vmax_vg > 0)
        vg_levels = [vg_levels, 0];  % ensure zero line if crossing
    end
    vg_levels = sort(vg_levels);

    fprintf('Meridional Wind (mid-lat): Max = %.2f, Min = %.2f [m/s]\n', vmax_vg, vmin_vg);

    %% --------------------------------------------------------------------
    %% 4. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4 2 18 14],'Visible','off');
    hold on;

    %% Filled contour: Geopotential Height
    contourf(xx, zz, gph_data', 'LineStyle','none');
    colormap(cmap_coolwarm(256));
    c1 = colorbar;
    c1.Label.String = 'Geopotential Height Perturbation (m)';
    c1.Label.FontSize = 20;
    caxis([vmin_gph vmax_gph]);

    %% Overlay: Contour lines for Meridional Wind
    [C, h] = contour(xx, zz, vg_data', vg_levels, 'LineColor','k', 'LineWidth',1.3);
    clabel(C, h, 'FontSize',14, 'Color','k', 'LabelSpacing',300);

    hold off;

    %% --------------------------------------------------------------------
    %% 5. Axes, labels, title
    %% --------------------------------------------------------------------
    xlabel('Longitude (degrees)', 'FontSize',20);
    ylabel('Height (km)', 'FontSize',20);
    set(gca, 'XTick', 0:30:360, 'YTick', 0:2:10);

    title_str = sprintf(['Geopotential Height (shaded) & Meridional Wind (contours)\n' ...
                         'Mid-Latitude Cross-Section; zonal wave # = %d, eMode # = %d'], ...
                         m0, n_mode);
    title(title_str, 'FontSize',20);

    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 6. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ...
        sprintf('combined_gph_meridional_xsec_eMode-%d_m0-%d.png', n_mode, m0));
    
    fprintf('Saving combined GPH + Meridional Wind plot to: %s\n', outFile);
    
    saveas(gcf, outFile);
    close(gcf);

end
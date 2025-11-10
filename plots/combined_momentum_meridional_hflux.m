%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_momentum_meridional_hflux.m

% DESCRIPTION: Overlays meridional heat flux <v'T'> as black contour lines
%              on a filled contour plot of momentum flux <v'u'>.

% INPUT:
%   vg     - 3-D meridional wind perturbation (m/s)
%   ug     - 3-D zonal wind perturbation (m/s)
%   temp   - 3-D temperature perturbation (K)
%   m0     - zonal wavenumber
%   n_mode - eigenmode index
%   fig_path - directory to save figure

% OUTPUT:
%   Saves: fig_path/combined_momentum_meridional_hflux_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_momentum_meridional_hflux(vg, ug, temp, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Compute zonal means
    %% --------------------------------------------------------------------
    zvu_calc = squeeze(mean(vg .* ug, 1));        % <v'u'>  (jj+1 × kk+1)
    zvt_calc = squeeze(mean(vg .* temp, 1));      % <v'T'>  (jj+1 × kk+1)

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4 2 18 14],'Visible','off');
    hold on;

    %% --------------------------------------------------------------------
    %% 3. Filled contour: Momentum flux <v'u'>
    %% --------------------------------------------------------------------
    contourf(zvu_calc', 'LineStyle','none');
    colormap(cmap_twilight(256));
    c1 = colorbar;
    c1.Label.String = '<v''u''> (m^2/s^2)';
    c1.Label.FontSize = 20;

    %% --------------------------------------------------------------------
    %% 4. Overlay: Contour lines for meridional heat flux <v'T'>
    %% --------------------------------------------------------------------
    [C, h] = contour(zvt_calc', 'LineColor','k','LineWidth',1.3);
    clabel(C, h, 'FontSize',14, 'Color','k', 'LabelSpacing',300);

    %% --------------------------------------------------------------------
    %% 5. Axes, title, fonts
    %% --------------------------------------------------------------------
    xlabel('Latitude','FontSize',20);
    ylabel('z-level','FontSize',20);

    title_str = sprintf(['Momentum Flux (shaded) & Meridional Heat Flux (contours)\n' ...
                         '(zonal wave # = %d, eMode # = %d)'], m0, n_mode);
    title(title_str,'FontSize',20);

    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 6. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ...
        sprintf('combined_momentum_meridional_hflux_eMode-%d_m0-%d.png', n_mode, m0));
    
    fprintf('\nSaving combined momentum + meridional heat flux plot to: %s\n', outFile);
    
    saveas(gcf, outFile);
    close(gcf);
end
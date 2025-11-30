%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_momentum_vertical_hflux.m

% DESCRIPTION: Overlays vertical heat flux <w'T'> as black contour lines
%              on a filled contour plot of momentum flux <v'u'>.

% INPUT:
%   vg      - 3-D meridional wind perturbation (m/s)
%   ug      - 3-D zonal wind perturbation (m/s)
%   wfield  - 3-D vertical velocity (m/s)
%   temp    - 3-D temperature perturbation (K)
%   m0      - zonal wavenumber
%   n_mode  - eigenmode index
%   fig_path- directory to save figure

% OUTPUT:
%   Saves: fig_path/combined_flux_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_momentum_vertical_hflux(vg, ug, wfield, temp, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Compute zonal means
    %% --------------------------------------------------------------------
    zvu_calc = squeeze(mean(vg .* ug, 1));        % <v'u'>  (jj+1 × kk+1)
    zwt_calc = squeeze(mean(wfield .* temp, 1));  % <w'T'>  (jj+1 × kk+1)

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
    c1.Label.String = '<v''u''> (m²/s²)';
    c1.Label.FontSize = 20;

    %% --------------------------------------------------------------------
    %% 4. Overlay: Contour lines for vertical heat flux <w'T'>
    %% --------------------------------------------------------------------
    [C, h] = contour(zwt_calc', 'LineColor','k','LineWidth',1.2);
    clabel(C, h, 'FontSize',14, 'Color','k','LabelSpacing',300);

    %% --------------------------------------------------------------------
    %% 5. Axes, title, fonts
    %% --------------------------------------------------------------------
    xlabel('Latitude','FontSize',20);
    ylabel('z-level','FontSize',20);

    title_str = sprintf(['Momentum Flux (shaded) & Vertical Heat Flux (contours)\n' ...
                         'zonal wave # = %d, eMode # = %d'], m0, n_mode);
    title(title_str,'FontSize',20);

    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 6. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ...
        sprintf('combined_flux_eMode-%d_m0-%d.png', n_mode, m0));
    fprintf('\nSaving combined momentum + vertical flux plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end
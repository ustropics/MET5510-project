%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_dpvdym_int.m

% DESCRIPTION: Plots the interior meridional gradient of background potential
%              vorticity (∂q̄/∂y) scaled by 10¹² s⁻¹ in the latitude-height
%              plane using filled contours. The figure is saved as a PNG
%              including wavenumber and mode number in title and filename.

% INPUT:
%   yy       - Latitude coordinates (degrees), size (jj+1)
%   zz       - Height coordinates (km), size (kk+1)
%   BPVy     - Meridional PV gradient field (jj+1 x kk+1 array, s⁻¹)
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

% OUTPUT:
%   Saves: fig_path/dpvdym-int_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dpvdym_int(yy, zz, BPVy, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Scale data and compute color limits
    %% --------------------------------------------------------------------
    data = BPVy * 1e12;                    % scale to 10⁻¹² s⁻¹
    vmin = min(data(:));
    vmax = max(data(:));
    
    fprintf('\nInterior meridional gradient of background potential vorticity:\n')
    fprintf('dPVdy (interior) - Max: %.3e, Min: %.3e [×10⁻¹² s⁻¹]\n', vmax, vmin);
    
    step = 0.2;
    vmin = floor(vmin / step) * step;
    vmax = ceil(vmax / step) * step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4,2,16,12],'Visible','off');
    contourf(yy, zz, data', 'LineStyle','none');
    hold on
    contour(yy, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off

    colormap(cmap_spectral(256));
    colorbar;
    caxis([vmin vmax]);

    xlabel('Latitude (degrees)');
    ylabel('Height (km)');
    
    title_str = ['d(PVbar)/dy (Interior, 10^-12 s^-1)', newline ...
    'zonal wave # = ', num2str(m0), ...
    ', eMode # = ', num2str(n_mode)];

    title(title_str);

    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ['dpvdym-int', ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving interior PV gradient plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end
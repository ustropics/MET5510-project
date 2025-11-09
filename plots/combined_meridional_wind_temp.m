%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_meridional_wind_temp.m

% DESCRIPTION: Plots meridional wind (filled contour) at a specified vertical
%              level across longitude and latitude, with temperature
%              overlaid as black contour lines (solid for positive,
%              dashed for negative, thin solid for zero). The figure is
%              saved as a PNG including hlevel, wavenumber, and mode number
%              in title and filename.

% INPUT:
%   xx       - Longitude coordinates (degrees)
%   yy       - Latitude coordinates (degrees)
%   vg       - Meridional wind field (ii+1 x jj+1 x kk+1 array, m/s)
%   temp     - Temperature field (same size as vg, K)
%   hlevel   - Vertical level index (e.g., 1, 25, 51)
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_meridional_wind_temp(xx, yy, vg, temp, hlevel, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4 2 16 12],'Visible','off')
    
    %% --------------------------------------------------------------------
    %% 2. Filled contour of meridional wind (background)
    %% --------------------------------------------------------------------
    contourf(xx, yy, squeeze(vg(:,:,hlevel))', 'LineStyle','none');
    colorbar;
    colormap(cmap_PRGn(256));
    hold on;
    
    %% --------------------------------------------------------------------
    %% 3. Overlay temperature contours: black, solid(+), dashed(-)
    %% --------------------------------------------------------------------
    temp_slice = squeeze(temp(:,:,hlevel))';
    levs = -0.4:0.1:0.4;                    % fine levels for small range
    zero_idx = find(levs == 0);
    
    % Positive: solid red
    if zero_idx < length(levs)
        contour(xx, yy, temp_slice, levs(zero_idx+1:end), 'r-', 'LineWidth', 2);
    end
    
    % Negative: dashed blue
    if zero_idx > 1
        contour(xx, yy, temp_slice, levs(1:zero_idx-1), 'b--', 'LineWidth', 2);
    end
    
    % Zero line: thin solid black (optional)
    contour(xx, yy, temp_slice, [0 0], 'k-', 'LineWidth', 0.8);

    %% --------------------------------------------------------------------
    %% 4. Axes, title, fonts
    %% --------------------------------------------------------------------
    hold off;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',0:30:360)
    
    title_str = ['Meridional Wind (shaded) with Temperature Differences (contoured)', newline ...
                 'hlevel = ', num2str(hlevel), ...
                 ', zonal wave # = ', num2str(m0), ...
                 ', eMode # = ', num2str(n_mode)];
    title(title_str);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    
    %% --------------------------------------------------------------------
    %% 5. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ['meridional_wind_temp', ...
                 '_hlevel-', num2str(hlevel), ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);

    fprintf('\nSaving meridional wind + temperature plot to: %s\n', outFile);
    
    saveas(gcf, outFile);
    close(gcf);
end
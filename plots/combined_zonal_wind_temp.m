%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_zonal_wind.m
%
% DESCRIPTION: Plots zonal wind (filled contour) at a specified vertical
%              level across longitude and latitude, with temperature
%              overlaid as black contour lines (solid for positive,
%              dashed for negative, thin solid for zero). The figure is
%              saved as a PNG including hlevel, wavenumber, and mode number
%              in title and filename.
%
% INPUT:
%   xx       - Longitude coordinates (degrees)
%   yy       - Latitude coordinates (degrees)
%   ug       - Zonal wind field (ii+1 x jj+1 x kk+1 array, m/s)
%   temp     - Temperature field (same size as ug, K)
%   hlevel   - Vertical level index (e.g., 1, 25, 51)
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_zonal_wind_temp(xx, yy, ug, temp, hlevel, m0, n_mode, fig_path)
    %% Create figure
    figure('units','inch','position',[4,2,16,12],'Visible','off')
    
    % Wind shading
    contourf(xx, yy, squeeze(ug(:,:,hlevel))', 'LineStyle','none');
    colorbar;
    colormap(cmap_RdYlBl(256));
    hold on;
    
    % --- Temperature contours: black, solid(+), dashed(-) ---
    temp_slice = squeeze(temp(:,:,hlevel))';
    levs = -0.6:0.1:0.6;                    % same fine levels as meridional plot
    zero_idx = find(levs == 0);
    
    % Positive: solid black
    if zero_idx < length(levs)
        contour(xx, yy, temp_slice, levs(zero_idx+1:end), 'k-', 'LineWidth', 2);
    end
    
    % Negative: dashed black
    if zero_idx > 1
        contour(xx, yy, temp_slice, levs(1:zero_idx-1), 'k--', 'LineWidth', 2);
    end
    
    % Zero line: thin solid black (optional — comment out to hide)
    contour(xx, yy, temp_slice, [0 0], 'k-', 'LineWidth', 0.8);
    % ---------------------------------------------------------
    
    hold off;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',0:30:360)
    title_str = ['Zonal Wind (hlevel = ', num2str(hlevel), ...
                 ', zonal wave # = ', num2str(m0), ...
                 ', eMode # = ', num2str(n_mode), ')'];
    title(title_str);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    
    %% Save
    outFile = fullfile(fig_path, ['zonal_wind_', ...
                 '_hlevel-', num2str(hlevel), ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);
    fprintf('Saving zonal wind + temperature plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end
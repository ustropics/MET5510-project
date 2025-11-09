%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_zonal_wind_temp.m

% DESCRIPTION: Plots zonal wind (filled contour) at a specified vertical
%              level across longitude and latitude, with temperature
%              overlaid as black contour lines (solid for positive,
%              dashed for negative, thin solid for zero).

% INPUT:
%   xx       - Longitude coordinates (degrees)
%   yy       - Latitude coordinates (degrees)
%   ug       - Zonal wind field (ii+1 x jj+1 x kk+1 array, m/s)
%   temp     - Temperature field (same size as ug, K)
%   hlevel   - Vertical level index (e.g., 1, 25, 51)
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

% OUTPUT:
%   Saves: fig_path/zonal_wind_temp_hlevel-<h>_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_zonal_wind_temp(xx, yy, ug, temp, hlevel, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Create figure
    %% --------------------------------------------------------------------

    figure('units','inch','position',[4,2,16,12],'Visible','off')
    
    %% --------------------------------------------------------------------
    %% 2. Filled contour of zonal wind (background)
    %% --------------------------------------------------------------------
    
    contourf(xx, yy, squeeze(ug(:,:,hlevel))', 'LineStyle','none');
    colorbar;
    colormap(cmap_PRGn(256));
    hold on;

    %% --------------------------------------------------------------------
    %% 3. Overlay temperature contours: black, solid(+), dashed(-)
    %% --------------------------------------------------------------------
   
    % Get temperature slice data and set levels
    temp_slice = squeeze(temp(:,:,hlevel))';
    levs = -1:0.1:1;  
    zero_idx = find(levs == 0);
    
    % Positive: solid black
    if zero_idx < length(levs)
        contour(xx, yy, temp_slice, levs(zero_idx+1:end), 'r-', 'LineWidth', 2);
    end
    
    % Negative: dashed black
    if zero_idx > 1
        contour(xx, yy, temp_slice, levs(1:zero_idx-1), 'b--', 'LineWidth', 2);
    end
    
    % Zero line: thin solid black (optional — comment out to hide)
    contour(xx, yy, temp_slice, [0 0], 'k-', 'LineWidth', 0.8);

    %% --------------------------------------------------------------------
    %% 4. Axes, title, fonts
    %% --------------------------------------------------------------------
   
    % Set x and y-label text as well as tick spacing for grid labels
    hold off;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',0:30:360)
    
    % Create title string with input variables and set it 
    title_str = ['Zonal Wind (shaded) with Temperature Differences (contoured)', newline ...
        'hlevel = ', num2str(hlevel), ...
        ', zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode)];
    
    title(title_str);

    % Set global font size    
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    
    %% --------------------------------------------------------------------
    %% 5. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['zonal_wind_temp', ...
                 '_hlevel-', num2str(hlevel), ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);

    fprintf('\nSaving zonal wind + temperature plot to: %s\n', outFile);
    
    % Comment out 'close(gcf)' if you want to see the plot and save it
    saveas(gcf, outFile);
    close(gcf);
end
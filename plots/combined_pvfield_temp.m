%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: combined_pvfield_temp.m

% DESCRIPTION: Plots potential vorticity (filled contour) at a specified vertical
%              level across longitude and latitude, with temperature
%              overlaid as colored contour lines (solid red for positive,
%              dashed blue for negative, thin solid black for zero). 

% INPUT:
%   xx       - Longitude coordinates (degrees)
%   yy       - Latitude coordinates (degrees)
%   pvfield  - Potential vorticity field (ii+1 x jj+1 x kk+1 array, s⁻¹)
%   temp     - Temperature field (same size as pvfield, K)
%   hlevel   - Vertical level index (e.g., 1 for surface)
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_pvfield_temp(xx, yy, pvfield, temp, hlevel, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4 2 16 12],'Visible','off')
    
    %% --------------------------------------------------------------------
    %% 2. Filled contour of potential vorticity (background)
    %% --------------------------------------------------------------------
    pv_slice = squeeze(pvfield(:,:,hlevel));
    contourf(xx, yy, pv_slice', 'LineStyle','none');
    
    % Symmetric color scaling (same as plot_pvfield.m)
    step = 0.2;
    vmin = floor(min(pv_slice(:))/step)*step;
    vmax = ceil (max(pv_slice(:))/step)*step;
    caxis([vmin vmax]);
    
    colorbar;
    colormap(cmap_coolwarm(256));  % or use your preferred colormap
    hold on;
    
    %% --------------------------------------------------------------------
    %% 3. Overlay temperature contours: red(+), blue(--), black zero
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
    
    % Zero line: thin solid black
    contour(xx, yy, temp_slice, [0 0], 'k-', 'LineWidth', 0.8);

    %% --------------------------------------------------------------------
    %% 4. Axes, title, fonts
    %% --------------------------------------------------------------------
    hold off;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',0:30:360)
    
    title_str = ['Potential Vorticity (shaded) with Temperature Differences (contoured)', newline ...
                 'hlevel = ', num2str(hlevel), ...
                 ', zonal wave # = ', num2str(m0), ...
                 ', eMode # = ', num2str(n_mode)];
    title(title_str);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    
    %% --------------------------------------------------------------------
    %% 5. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ['pvfield_temp', ...
                 '_hlevel-', num2str(hlevel), ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);

    fprintf('\nSaving potential vorticity + temperature plot to: %s\n', outFile);
    
    saveas(gcf, outFile);
    close(gcf);
end
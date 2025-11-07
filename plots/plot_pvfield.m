%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_pvfield.m

% DESCRIPTION: Plots the potential vorticity contour at the surface (k=1)
% across longitude and latitude, and saves it as an image.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - pvfield: Potential vorticity field (ii+1 x jj+1 x kk+1 array, s^-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_pvfield(xx, yy, pvfield, m0, n_mode, fig_path)

    %% create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(pvfield(:,:,1))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)

    % set the title
    title_str = ['Potential Vorticity (', ...
    'zonal wave # = ', num2str(m0), ...
    ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% save figure
    outFile = fullfile(fig_path, ['pvfield', '_eMode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving temperature plot to: %s\n', outFile)
    saveas(gcf, outFile);
    close(gcf);
    
end
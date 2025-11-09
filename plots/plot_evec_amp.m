%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_evec_amp.m

% DESCRIPTION: Plots eigenvector amplitude in the latitude-height plane
%              using filled contours. Includes growth rate and frequency
%              in title. The figure is saved as a PNG with wavenumber and
%              mode number in filename.

% INPUT:
%   yy          - Latitude coordinates (degrees), size (jj+1)
%   zz          - Height coordinates (km), size (kk+1)
%   eVec_amp    - Eigenvector amplitude (jj+1 x kk+1 array)
%   m0          - Zonal wavenumber
%   n_mode      - Mode number
%   growth_rate - Real part of eigenvalue (growth rate)
%   omega       - Imaginary part of eigenvalue (frequency)
%   fig_path    - Directory path for saving figure

% OUTPUT:
%   Saves: fig_path/evec-amp_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    % Get 2d data for our eigenvector amplitude contour
    data = eVec_amp; % get our data to plot

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nEigenvector amplitude contour:\n')
    fprintf('Maximum Value: %.2e and Minimum Value: %.0f\n', vmax, vmin)
    
    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------

    % Create the figure and set it's size [left, bottom, width, height]
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(yy, zz, data', 'linestyle', 'none');
    hold on
    contour(yy, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off

    colormap(flipud(cmap_spectral(256))); % Set our colormap choice
    colorbar;
    hold on

    % Set x and y-label text as well as tick spacing for grid labels
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')
    set(gca, 'fontsize', 22, 'color', 'w')

    % Create title string with input variables and set it
    title_str = ['Zonal wave # = ', num2str(m0), ...
        ', eMode # =', num2str(n_mode), ...
        ', real eVal = ', num2str(growth_rate), ...
        ', imag eVal = ', num2str(omega)];

    title(title_str);
    
    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);


    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['evac-amp', ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving eigenvector amplitude plot to: %s\n', outFile)

    % Comment out 'close(gcf)' if you want to see the plot and save it    
    saveas(gcf, outFile);
    close(gcf);

end
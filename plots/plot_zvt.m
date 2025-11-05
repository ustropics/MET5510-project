function plot_zvt(vg, temp, m0, n_mode, fig_path)

    %% calculate zvt
    zvt_calc = squeeze(mean(vg.*temp,1));

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(zvt_calc', 'LineStyle', 'none');
    colorbar;
    xlabel('Latitude')
    ylabel('z-level')
    % set(gca, 'xtick', 0:30:360)

    title_str = ['Mean Meridional Temperature', ' (zonal wave # = ', ...
     num2str(m0), ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

        
    %% Save figure
    outFile = fullfile(fig_path, ['mean_meridional_temp', '_nmode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    fprintf('Saving zonal wind plot to: %s\n', outFile); % Debug output
    saveas(gcf, outFile);
    close(gcf);
    
end

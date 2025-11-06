function plot_zwt(wfield, temp, m0, n_mode, fig_path)

    %% calculate zvt
    zwt_calc = squeeze(mean(wfield.*temp,1));

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(zwt_calc, 'LineStyle', 'none');
    colorbar;
    xlabel('Latitude')
    ylabel('z-level')
    % set(gca, 'xtick', 0:30:360)

    title_str = ['Mean Vertical Temperature Advection', ' (zonal wave # = ', ...
     num2str(m0), ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

        
    %% Save figure
    outFile = fullfile(fig_path, ['mean_vertical_Tadv', '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    fprintf('Saving zonal wind plot to: %s\n', outFile); % Debug output
    saveas(gcf, outFile);
    close(gcf);
    
end

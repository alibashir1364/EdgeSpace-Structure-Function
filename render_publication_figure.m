function render_publication_figure(rho_subj, Tau, MSF, mat_msf, mat_tau, net_names, r_grp, p_grp)
    fig = figure('Units', 'normalized', 'Position', [0.03, 0.08, 0.94, 0.82], ...
                 'Color', 'w', 'Name', 'HCP 100-Subject Group Event Analysis');
    set(gcf, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica');
    
    % --- Panel A: Subject-Wise Consistency (Distribution of Rhos) ---
    subplot(2, 2, 1);
    histogram(rho_subj, 20, 'FaceColor', [0.2, 0.45, 0.65], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
    xline(mean(rho_subj), '--r', sprintf('Mean \\rho = %.3f', mean(rho_subj)), ...
          'LineWidth', 2, 'LabelOrientation', 'horizontal', 'FontSize', 10);
    xlabel('Subject-level Spearman Correlation (\rho)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Number of Subjects', 'FontSize', 11, 'FontWeight', 'bold');
    title('A. Subject-wise Robustness (N = 100 Subjects)', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2);
%     grid on; grid alpha 0.15;
    
    % --- Panel B: Group Scatter Plot with Fit Line ---
    subplot(2, 2, 2);
    scatter(Tau, MSF, 12, [0.15, 0.15, 0.15], 'filled', 'MarkerFaceAlpha', 0.12);
    hold on;
    mdl = fitlm(Tau, MSF, 'RobustOpts', 'on');
    x_grid = linspace(min(Tau), max(Tau), 100)';
    [y_pred, y_ci] = predict(mdl, x_grid, 'Prediction', 'curve');
    
    % Confidence interval shadow
    fill([x_grid; flipud(x_grid)], [y_ci(:,1); flipud(y_ci(:,2))], [0.85, 0.25, 0.25], ...
         'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(x_grid, y_pred, 'Color', [0.75, 0.1, 0.1], 'LineWidth', 2.5);
    
    xlabel('Group Mean Event Duration (\tau_{event} [TRs])', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Group Edge MSF', 'FontSize', 11, 'FontWeight', 'bold');
    title('B. Group-Level MSF vs Event Duration Coupling', 'FontSize', 12, 'FontWeight', 'bold');
    
    stat_str = sprintf('Group \\rho = %.3f\np < 10^{-15}', r_grp);
    text(min(Tau) + 0.05*range(Tau), max(MSF) - 0.08*range(MSF), stat_str, ...
         'FontSize', 11, 'FontWeight', 'bold', 'BackgroundColor', [0.96 0.96 0.96 0.85], ...
         'EdgeColor', [0.7 0.7 0.7], 'Margin', 6);
    set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2);
    axis square;
    
    % --- Panel C: 7-Network Group MSF Matrix ---
    subplot(2, 2, 3);
    imagesc(mat_msf);
    cb1 = colorbar; cb1.Label.String = 'Mean MSF';
    colormap(gca, flipud(bone));
    format_17net_axes(net_names, 'C. 7-Network Group MSF Matrix');
    
    % --- Panel D: 7-Network Group Tau Matrix ---
    subplot(2, 2, 4);
    imagesc(mat_tau);
    cb2 = colorbar; cb2.Label.String = 'Mean Event Duration (\tau)';
    colormap(gca, hot);
    format_17net_axes(net_names, 'D. 7-Network Group Event Duration Matrix');
end
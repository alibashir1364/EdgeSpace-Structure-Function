% -------------------------------------------------------------------
% PANEL A: Horizontal raincloud (PUBLICATION-style)
% -------------------------------------------------------------------
rng(42);                                

axA = axes('Position',[0.06 0.16 0.20 0.72]); hold(axA,'on');


col_violin = [0.36 0.55 0.74];   % density fill
col_edge   = [0.15 0.15 0.15];   % edges
col_emp    = [0.10 0.20 0.35];   % points
col_med    = [0.85 0.33 0.10];   % median highlight


[f_k, x_k] = ksdensity(rho_s, 'NumPoints', 300);
f_k = f_k / max(f_k);                    
y_base = 0.25;                            
y_top  = 0.90;                           
f_scaled = y_base + f_k*(y_top - y_base);

patch([x_k, fliplr(x_k)], [f_scaled, y_base*ones(size(f_scaled))], ...
      col_violin, 'FaceAlpha', 0.40, ...
      'EdgeColor', col_edge, 'LineWidth', 0.6);


hb = boxplot(rho_s, 'Orientation','horizontal', ...
             'Positions', 0.12, 'Widths', 0.07, ...
             'Symbol','', 'Colors', col_edge, 'Whisker', 1.5);
set(findobj(hb,'Tag','Box'),    'LineWidth', 0.8);
set(findobj(hb,'Tag','Whisker'),'LineWidth', 0.8);
set(findobj(hb,'Tag','Median'), 'LineWidth', 1.4, 'Color', col_med);


y_j = -0.06 - 0.10*rand(numel(rho_s),1);
scatter(rho_s, y_j, 14, col_emp, 'filled', ...
        'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor','none');


xline(0, '--', 'Color',[0.55 0.55 0.55], 'LineWidth', 0.7);


xlabel('Subject-level coupling  \rho_s', ...
       'FontName','Arial','FontSize', 8);
title('a | Subject distribution', ...
      'FontName','Arial','FontSize', 9, 'FontWeight','bold');

xlim([-0.55 -0.15]); ylim([-0.22 1.0]);
set(axA, 'YTick',[], 'TickDir','out', 'Box','off', ...
    'FontName','Arial','FontSize', 8, 'LineWidth', 0.75, ...
    'XColor',[0.15 0.15 0.15], 'YColor',[0.15 0.15 0.15]);


txt = sprintf(['N = %d\nMean = %.3f\nSD = %.3f\n100%% negative\np_{Wilcoxon} < 0.001'], ...
              numel(rho_s), mean_rho, sd_rho);
text(-0.16, 0.92, txt, 'FontName','Arial','FontSize', 7, ...
     'HorizontalAlignment','right', 'VerticalAlignment','top', ...
     'BackgroundColor',[1 1 1 0.9], ...
     'EdgeColor',[0.75 0.75 0.75], 'LineWidth',0.5, 'Margin', 3);
%% panel c mixed-effect
%% ========================================================================
%  PANEL C — VARIANCE DECOMPOSITION (Nature-style)
%  Input: sigma_subj, sigma_res, icc, beta, tval, pval
%  Output: publication-ready figure (PNG + SVG + EPS)
% ========================================================================

if ~exist('sigma_subj','var'),  sigma_subj = 0.03207;  end
if ~exist('sigma_res','var'),   sigma_res  = 0.07205;  end
if ~exist('icc','var'),         icc        = 0.1654;   end
if ~exist('beta','var'),        beta       = -0.006956; end
if ~exist('tval','var'),        tval       = -176.99;  end
if ~exist('pval','var'),        pval       = 0;        end

% --- Compute variances & percentages ---
var_sub = sigma_subj^2;
var_res = sigma_res^2;
var_tot = var_sub + var_res;

pct_sub = 100 * var_sub / var_tot;   % 16.54
pct_res = 100 * var_res / var_tot;   % 83.46

% --- PUBLICATION/NPG-inspired colors ---
col_sub = [0.23, 0.33, 0.53];   % Navy blue
col_res = [0.85, 0.85, 0.85];   % Light grey
col_txt = [0.15, 0.15, 0.15];

% ========================================================================
%  FIGURE
% ========================================================================
fig = figure('Units','normalized', 'Position', [0.30, 0.30, 0.22, 0.45], ...
             'Color', 'w', 'Name', 'PanelC_VarianceDecomposition');
set(gcf, 'DefaultAxesFontName', 'Helvetica', ...
         'DefaultTextFontName', 'Helvetica');

ax = axes('Position', [0.20, 0.18, 0.72, 0.72]);
hold(ax, 'on');

% --- Draw bars -------------------------------------------------------
b = bar(ax, [1, 2], [pct_sub, pct_res], 0.55, ...
        'FaceColor', 'flat', ...
        'EdgeColor', [0.15, 0.15, 0.15], ...
        'LineWidth', 1.2);

b.CData(1, :) = col_sub;
b.CData(2, :) = col_res;

% --- Add percentage labels on top of bars ----------------------------
text(1, pct_sub + 2, sprintf('%.1f%%', pct_sub), ...
     'FontSize', 11, 'FontWeight', 'bold', 'Color', col_sub, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(2, pct_res + 2, sprintf('%.1f%%', pct_res), ...
     'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.35, 0.35, 0.35], ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% --- Add ICC annotation with arrow to between-subject bar ------------
annotation_y = 45;
text(1, annotation_y, sprintf('ICC = %.2f', icc), ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', col_sub, ...
     'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
     'EdgeColor', col_sub, 'LineWidth', 1, 'Margin', 4);

% Dashed connector line from text to bar
plot([1, 1], [annotation_y - 3, pct_sub + 2.5], '--', ...
     'Color', col_sub, 'LineWidth', 1);

% --- Statistical annotation box (top-left) ---------------------------
stat_str = sprintf(['\\beta = %.5f\n' ...
                    't(%d) = %.1f\n' ...
                    'p < 10^{-300}'], ...
                    beta, 199998, tval);

text(1.55, 68, stat_str, ...
     'FontSize', 10, 'FontName', 'Helvetica', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'BackgroundColor', [1 1 1 0.95], 'EdgeColor', [0.6 0.6 0.6], ...
     'LineWidth', 0.8, 'Margin', 6);

% --- Axis labels and limits ------------------------------------------
set(ax, 'XTick', [1, 2], ...
        'XTickLabel', {'Between-subject', 'Residual'}, ...
        'FontSize', 11, 'TickDir', 'out', ...
        'Box', 'off', 'LineWidth', 1.2, ...
        'XLim', [0.4, 2.6], ...
        'YLim', [0, 100]);

ylabel('Variance explained (%)', 'FontSize', 11.5, 'FontWeight', 'bold');

title('C | Variance decomposition', 'FontSize', 12, 'FontWeight', 'bold');

% --- Turn off grid & make it clean -----------------------------------
grid(ax, 'off');

% --- Optional: light horizontal guide at 50% -------------------------
plot([0.4, 2.6], [50, 50], ':', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.7);

% ========================================================================
%  EXPORT
% ========================================================================
saveas(fig, 'PanelC_VarianceDecomposition.png');
saveas(fig, 'PanelC_VarianceDecomposition.svg');
try
    exportgraphics(fig, 'PanelC_VarianceDecomposition.eps', 'ContentType', 'vector');
catch
    print(fig, 'PanelC_VarianceDecomposition.eps', '-depsc2', '-r300');
end

fprintf('Panel C saved: PNG, SVG, EPS\n');
fprintf('  Between-subject variance: %.2f%% (ICC = %.4f)\n', pct_sub, icc);
fprintf('  Residual variance:        %.2f%%\n', pct_res);
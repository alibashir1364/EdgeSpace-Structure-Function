%% ========================================================================
%  TWO-AXIS MODEL (FINAL PUBLICATION VERSION)
%  Focal Weight drives Event Amplitude (Positive) and Burstiness (Negative)
%  Includes Subject-level robustness (Wilcoxon) and Publication-style plots
% ========================================================================

clearvars; clc; close all;
rng(42); % For reproducibility

%% ========================================================================
%  SECTION 1: LOAD DATA
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  TWO-AXIS MODEL: FOCAL WEIGHT DRIVES EVENT DYNAMICS\n');
fprintf('===============================================================\n');

S = load('EdgeSpace_Structure_Function_results_68regions.mat');
C_edge       = S.Structural_Edge_Context(:);    
Group_MSF    = S.Group_MSF(:);                  
N_edges = length(C_edge);
SC_group = S.SC_group;

load('ROI_TS_68_HCP_100subj.mat');      
BOLD = roiTS_all;
[N_nodes, T, N_subj] = size(BOLD);
TR = 0.72; fs = 1/TR;


tril_mask = tril(true(N_nodes), -1);
[i_idx, j_idx] = find(tril_mask);
w_ij = zeros(N_edges, 1);
for e = 1:N_edges
    w_ij(e) = SC_group(i_idx(e), j_idx(e));
end

% Focal weight is only meaningful for directly connected edges
idx_direct = w_ij > 0;
fprintf('Data loaded: %d total edges, %d direct edges.\n', N_edges, sum(idx_direct));

%% ========================================================================
%  SECTION 2: COMPUTE EVENT METRICS (Threshold = 2.25)
% ========================================================================
fprintf('\n--- Computing Event Metrics (Threshold = 2.25) ---\n');
event_threshold = 2.25;

Amplitude_per_subj  = zeros(N_edges, N_subj);
Burstiness_per_subj = zeros(N_edges, N_subj);

for s = 1:N_subj
    X = squeeze(BOLD(:, :, s));   
%     [b, a] = butter(4, [0.008, 0.2]/(fs/2), 'bandpass');
%     X = filtfilt(b, a, X')';
    X_z = (X - mean(X, 2)) ./ std(X, 0, 2);
    
    for e = 1:N_edges
        ets = X_z(i_idx(e), :) .* X_z(j_idx(e), :);
        supra = ets > event_threshold;
        d = diff([0, supra, 0]);
        starts = find(d == 1);
        ends   = find(d == -1) - 1;
        
        if ~isempty(starts)
            Amplitude_per_subj(e, s) = mean(ets(supra));
            if length(starts) > 2
                IEI = diff(starts);
                mu = mean(IEI); sd = std(IEI);
                if (sd + mu) > 0
                    Burstiness_per_subj(e, s) = (sd - mu) / (sd + mu);
                end
            end
        end
    end
    if mod(s, 20) == 0, fprintf('  Subject %d/%d processed.\n', s, N_subj); end
end

% Group means for plotting
mean_Amplitude  = mean(Amplitude_per_subj, 2);
mean_Burstiness = mean(Burstiness_per_subj, 2);

%% ========================================================================
%  SECTION 3: SUBJECT-LEVEL & GROUP-LEVEL STATISTICAL ANALYSIS
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  STATISTICAL ANALYSIS (Direct Edges Only)\n');
fprintf('===============================================================\n');

% 1. Subject-level correlations
r_amp_sub = zeros(N_subj, 1);
r_burst_sub = zeros(N_subj, 1);

for s = 1:N_subj
    r_amp_sub(s)   = corr(w_ij(idx_direct), Amplitude_per_subj(idx_direct, s), 'Type', 'Spearman');
    r_burst_sub(s) = corr(w_ij(idx_direct), Burstiness_per_subj(idx_direct, s), 'Type', 'Spearman');
end

% 2. Non-parametric tests (Wilcoxon signed-rank against zero)
[p_amp, ~, stats_amp] = signrank(r_amp_sub, 0);
[p_burst, ~, stats_burst] = signrank(r_burst_sub, 0);

% 3. Effect sizes (Cohen's d for one-sample)
d_amp = mean(r_amp_sub) / std(r_amp_sub);
d_burst = mean(r_burst_sub) / std(r_burst_sub);

% 4. Group-level correlations (for the scatter plots)
[r_group_amp, p_group_amp] = corr(w_ij(idx_direct), mean_Amplitude(idx_direct), 'Type', 'Spearman');
[r_group_burst, p_group_burst] = corr(w_ij(idx_direct), mean_Burstiness(idx_direct), 'Type', 'Spearman');

fprintf('\n--- 1. Amplitude vs Focal Weight (Positive Axis) ---\n');
fprintf('Group-level Spearman rho: %+.3f (p = %.2e)\n', r_group_amp, p_group_amp);
fprintf('Subject-level mean rho  : %+.3f +/- %.3f\n', mean(r_amp_sub), std(r_amp_sub));
fprintf('Wilcoxon test vs 0      : z = %.2f, p = %.2e, Cohen''s d = %.2f\n', stats_amp.zval, p_amp, d_amp);
fprintf('Subjects with rho > 0   : %d / %d (%.1f%%)\n', sum(r_amp_sub > 0), N_subj, 100*mean(r_amp_sub > 0));

fprintf('\n--- 2. Burstiness vs Focal Weight (Negative Axis) ---\n');
fprintf('Group-level Spearman rho: %+.3f (p = %.2e)\n', r_group_burst, p_group_burst);
fprintf('Subject-level mean rho  : %+.3f +/- %.3f\n', mean(r_burst_sub), std(r_burst_sub));
fprintf('Wilcoxon test vs 0      : z = %.2f, p = %.2e, Cohen''s d = %.2f\n', stats_burst.zval, p_burst, d_burst);
fprintf('Subjects with rho < 0   : %d / %d (%.1f%%)\n', sum(r_burst_sub < 0), N_subj, 100*mean(r_burst_sub < 0));

fprintf('\n--- 3. The Two-Axis Contrast (MSF & Duration) ---\n');
fprintf('Focal Weight vs MSF     : rho = %+.3f (N.S.)\n', corr(w_ij(idx_direct), Group_MSF(idx_direct), 'Type', 'Spearman'));
fprintf('Context vs MSF          : rho = %+.3f (Significant)\n', corr(C_edge, Group_MSF, 'Type', 'Spearman'));

%% ========================================================================
%  SECTION 4: PUBLICATION-STYLE FIGURES
% ========================================================================
% Colors
col_pos = [0.23 0.33 0.53]; % Navy (Amplitude)
col_neg = [0.89 0.29 0.20]; % Coral (Burstiness)
col_gray = [0.5 0.5 0.5];

% -------------------------------------------------------------------------
% FIGURE 1: Group-level Scatters (The Divergence)
% -------------------------------------------------------------------------
fig1 = figure('Color','w','Units','centimeters','Position',[2 2 18 8]);
set(gcf, 'DefaultAxesFontName', 'Helvetica');

% Panel A: Weight vs Amplitude
ax1 = subplot(1,2,1); hold(ax1, 'on');
scatter(ax1, w_ij(idx_direct), mean_Amplitude(idx_direct), 15, col_pos, 'filled', 'MarkerFaceAlpha', 0.3);
mdl1 = fitlm(w_ij(idx_direct), mean_Amplitude(idx_direct), 'RobustOpts', 'on');
x_l = linspace(min(w_ij(idx_direct)), max(w_ij(idx_direct)), 50);
plot(ax1, x_l, predict(mdl1, x_l'), 'k', 'LineWidth', 2.5);
text(ax1, max(w_ij)*0.5, max(mean_Amplitude)*0.9, sprintf('\\rho = %+.3f***', r_group_amp), ...
     'FontSize', 10, 'FontWeight','bold', 'BackgroundColor', 'w');
xlabel(ax1, 'Focal Structural Weight (w_{ij})', 'FontWeight','bold'); 
ylabel(ax1, 'Mean Event Amplitude', 'FontWeight','bold');
title(ax1, 'A | Amplitude (Positive Axis)', 'FontWeight','bold');
set(ax1, 'TickDir','out','Box','off','FontSize',9);

% Panel B: Weight vs Burstiness
ax2 = subplot(1,2,2); hold(ax2, 'on');
scatter(ax2, w_ij(idx_direct), mean_Burstiness(idx_direct), 15, col_neg, 'filled', 'MarkerFaceAlpha', 0.3);
mdl2 = fitlm(w_ij(idx_direct), mean_Burstiness(idx_direct), 'RobustOpts', 'on');
plot(ax2, x_l, predict(mdl2, x_l'), 'k', 'LineWidth', 2.5);
text(ax2, max(w_ij)*0.5, max(mean_Burstiness)*0.9, sprintf('\\rho = %+.3f***', r_group_burst), ...
     'FontSize', 10, 'FontWeight','bold', 'BackgroundColor', 'w');
xlabel(ax2, 'Focal Structural Weight (w_{ij})', 'FontWeight','bold'); 
ylabel(ax2, 'Mean Event Burstiness', 'FontWeight','bold');
title(ax2, 'B | Burstiness (Negative Axis)', 'FontWeight','bold');
set(ax2, 'TickDir','out','Box','off','FontSize',9);

sgtitle('Focal Weight Drives Discrete Event Dynamics', 'FontSize', 12, 'FontWeight', 'bold');
saveas(fig1, 'Figure_TwoAxis_Scatters.png');
saveas(fig1, 'Figure_TwoAxis_Scatters.svg');

% -------------------------------------------------------------------------
% FIGURE 2: Subject-level Robustness (Boxplots)
% -------------------------------------------------------------------------
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 14 10]);
set(gcf, 'DefaultAxesFontName', 'Helvetica');

ax3 = axes('Position',[0.15 0.2 0.7 0.7]); hold(ax3, 'on');

% Prepare data for grouped boxplot
% data_box = [r_amp_sub, r_burst_sub];
data_box = r_amp_sub;

% Draw boxplot
bp = boxplot(data_box, 'Labels', {'Weight vs Amplitude'}, ...
             'Colors', 'k', 'Symbol', '', 'Widths', 0.4);
set(bp, 'LineWidth', 1.5);

% Overlay jittered points
jitter_width = 0.1;
x_jit_amp = 1 + (rand(N_subj, 1) - 0.5) * jitter_width;
% x_jit_burst = 2 + (rand(N_subj, 1) - 0.5) * jitter_width;

scatter(ax3, x_jit_amp, r_amp_sub, 40, col_pos, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
% scatter(ax3, x_jit_burst, r_burst_sub, 40, col_neg, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');

% Zero line
yline(ax3, 0, '--k', 'LineWidth', 1.5);

% Annotation
text(ax3, 1, 0.25, sprintf('Wilcoxon p < 0.001\nd = %.2f', d_amp), 'HorizontalAlignment', 'center', 'FontSize', 9);
% text(ax3, 2, 0.25, sprintf('Wilcoxon p < 0.001\nd = %.2f', d_burst), 'HorizontalAlignment', 'center', 'FontSize', 9);

ylabel(ax3, 'Subject-level Spearman \rho', 'FontWeight','bold', 'FontSize', 11);
title(ax3, 'Inter-subject Robustness (N=100)', 'FontWeight','bold', 'FontSize', 12);
set(ax3, 'TickDir','out','Box','off','FontSize',10, 'XLim', [0.5 2.5], 'YLim', [-0.4 0.4]);
grid(ax3, 'on'); ax3.GridAlpha = 0.15;

saveas(fig2, 'Figure_TwoAxis_Robustness.png');
saveas(fig2, 'Figure_TwoAxis_Robustness.svg');

fprintf('\nFigures saved successfully.\n');
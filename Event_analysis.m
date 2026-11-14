%% ========================================================================
%  BRIDGING SPECTRAL AND EVENT-BASED VIEWS OF EDGE TIMESCALES
%  Following Faskowitz et al. (2026) for event extraction
%  + Our MSF/Context framework for cross-validation
% ========================================================================

clearvars; clc; close all;

%% ========================================================================
%  SECTION 1: LOAD DATA
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  BRIDGING SPECTRAL AND EVENT-BASED EDGE TIMESCALES\n');
fprintf('===============================================================\n');

% --- Load structural edge-space results
S = load('EdgeSpace_Structure_Function_results_68regions.mat');
C_edge       = S.Structural_Edge_Context(:);    
Group_MSF    = S.Group_MSF(:);                  
mask_direct  = logical(S.Direct_mask(:));
mask_nondir  = logical(S.NonDirect_mask(:));

N_edges = length(C_edge);
fprintf('Edges loaded: N = %d\n', N_edges);

% --- Load BOLD time series
% Adjust the file name if needed
load('ROI_TS_68_HCP_100subj.mat');    
BOLD = roiTS_all;
[N_nodes, T, N_subj] = size(BOLD);
TR = 0.72;   % Repetition time (seconds)
fs = 1 / TR;

fprintf('BOLD: %d nodes x %d time points x %d subjects\n', N_nodes, T, N_subj);

% --- Reconstruct edge indices (lower triangle) ---
tril_mask = tril(true(N_nodes), -1);
[i_idx, j_idx] = find(tril_mask);

%% ========================================================================
%  SECTION 2: COMPUTE EDGE TIME SERIES AND EVENT DURATION
%  Following Faskowitz et al. (2026): threshold = 2.25
% ========================================================================
fprintf('\n--- Computing edge time series and event durations ---\n');

event_threshold = 2.25;   % Following Faskowitz et al. (2026)
% event_threshold = 1.5;
min_duration    = 1;      % Minimum frames to count as event
max_duration    = 100;    % Cap to avoid infinite loops

% Pre-allocate: mean event duration per edge per subject
Duration_per_subj = zeros(N_edges, N_subj);
EventCount_per_subj = zeros(N_edges, N_subj);

for s = 1:N_subj
    % --- Step 1: Z-score each node time series ---
    X = squeeze(BOLD(:, :, s));   % N_nodes x T
    % 1. Bandpass 0.008 - 0.2 Hz
    [b, a] = butter(4, [0.008, 0.2]/(fs/2), 'bandpass');
    X = filtfilt(b, a, X')';
    X = detrend(X')';
    X_z = (X - mean(X, 2)) ./ std(X, 0, 2);
    
    % --- Step 2: Compute edge time series (element-wise product) ---
    % For each edge (i,j): ETS(t) = X_z(i,t) * X_z(j,t)
    for e = 1:N_edges
        ets = X_z(i_idx(e), :) .* X_z(j_idx(e), :);
        
        % --- Step 3: Threshold ---
        supra = ets > event_threshold;
        
        % --- Step 4: Identify contiguous supra-threshold runs ---
        % Find transitions: 0->1 (start) and 1->0 (end)
        d = diff([0, supra, 0]);
        starts = find(d == 1);
        ends   = find(d == -1) - 1;
        durations = ends - starts + 1;
        
        % --- Step 5: Compute mean duration ---
        valid_dur = durations(durations >= min_duration & durations <= max_duration);
        
        if isempty(valid_dur)
            Duration_per_subj(e, s) = 0;
            EventCount_per_subj(e, s) = 0;
        else
            Duration_per_subj(e, s) = mean(valid_dur) * TR;  % Convert to seconds
            EventCount_per_subj(e, s) = length(valid_dur);
        end
    end
    
    if mod(s, 10) == 0
        fprintf('  Subject %d/%d processed\n', s, N_subj);
    end
end

% --- Group-level mean duration (in seconds) ---
Duration_group_mean = mean(Duration_per_subj, 2);
EventCount_group_mean = mean(EventCount_per_subj, 2);

fprintf('Event duration computed: mean = %.2f s (range: %.2f - %.2f s)\n', ...
        mean(Duration_group_mean), min(Duration_group_mean), max(Duration_group_mean));

%% ========================================================================
%  SECTION 3: CORRELATION ANALYSES
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  CORRELATION ANALYSES\n');
fprintf('===============================================================\n');

% --- Filter edges with valid duration (> 0) ---
valid_edges = Duration_group_mean > 0;
fprintf('Valid edges: %d/%d\n', sum(valid_edges), N_edges);

C_valid = C_edge(valid_edges);
MSF_valid = Group_MSF(valid_edges);
Dur_valid = Duration_group_mean(valid_edges);

% --- (A) MSF vs Event Duration ---
[rho_A, p_A] = corr(MSF_valid, Dur_valid, 'Type', 'Spearman');
fprintf('\n(A) MSF vs Event Duration:\n');
fprintf('    rho = %+.3f, p = %.2e, R2 = %.3f\n', rho_A, p_A, rho_A^2);

% --- (B) Structural Context vs Event Duration ---
[rho_B, p_B] = corr(C_valid, Dur_valid, 'Type', 'Spearman');
fprintf('\n(B) Structural Context vs Event Duration:\n');
fprintf('    rho = %+.3f, p = %.2e, R2 = %.3f\n', rho_B, p_B, rho_B^2);

% --- (C) Event Count vs Event Duration (control) ---
Count_valid = EventCount_group_mean(valid_edges);
[rho_C, p_C] = corr(Count_valid, Dur_valid, 'Type', 'Spearman');
fprintf('\n(C) Event Count vs Event Duration (control):\n');
fprintf('    rho = %+.3f, p = %.2e, R2 = %.3f\n', rho_C, p_C, rho_C^2);

% --- (D) Partial correlation: Context vs Duration | MSF ---
[rho_partial, p_partial] = partialcorr(C_valid, Dur_valid, MSF_valid, ...
                                        'Type', 'Spearman');
fprintf('\n(D) Partial correlation: Context vs Duration | MSF:\n');
fprintf('    rho_partial = %+.3f, p = %.2e\n', rho_partial, p_partial);

% --- (E) Direct vs Non-direct edges ---
direct_valid    = mask_direct(valid_edges);
nondirect_valid = mask_nondir(valid_edges);

if sum(direct_valid) > 0 && sum(nondirect_valid) > 0
    [rho_dir, p_dir] = corr(C_valid(direct_valid), Dur_valid(direct_valid), ...
                             'Type', 'Spearman');
    [rho_non, p_non] = corr(C_valid(nondirect_valid), Dur_valid(nondirect_valid), ...
                             'Type', 'Spearman');
    fprintf('\n(E) Direct vs Non-direct:\n');
    fprintf('    Direct:     rho = %+.3f, p = %.2e (N = %d)\n', ...
            rho_dir, p_dir, sum(direct_valid));
    fprintf('    Non-direct: rho = %+.3f, p = %.2e (N = %d)\n', ...
            rho_non, p_non, sum(nondirect_valid));
end

%% ========================================================================
%  SECTION 4: BOOTSTRAP VALIDATION FOR CORRELATIONS
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  BOOTSTRAP VALIDATION\n');
fprintf('===============================================================\n');

N_boot = 5000;
boot_A = zeros(N_boot, 1);
boot_B = zeros(N_boot, 1);
rng(42);

n_valid = sum(valid_edges);
for b = 1:N_boot
    idx = randsample(n_valid, n_valid, true);
    boot_A(b) = corr(MSF_valid(idx), Dur_valid(idx), 'Type', 'Spearman');
    boot_B(b) = corr(C_valid(idx),   Dur_valid(idx), 'Type', 'Spearman');
end

ci_A = prctile(boot_A, [2.5, 97.5]);
ci_B = prctile(boot_B, [2.5, 97.5]);

fprintf('MSF vs Duration:         95%% CI = [%+.3f, %+.3f]\n', ci_A(1), ci_A(2));
fprintf('Context vs Duration:     95%% CI = [%+.3f, %+.3f]\n', ci_B(1), ci_B(2));

%% ========================================================================
%  SECTION 5: SPATIAL MAPPING (NODE-LEVEL)
% ========================================================================
fprintf('\n--- Computing node-level aggregation ---\n');

% Aggregate edge-level metrics to node level
node_duration = zeros(N_nodes, 1);
node_MSF      = zeros(N_nodes, 1);
node_context  = zeros(N_nodes, 1);

for n = 1:N_nodes
    % Edges emanating from node n
    edge_from_n = (i_idx == n) | (j_idx == n);
    
    node_duration(n) = mean(Duration_group_mean(edge_from_n & valid_edges), 'omitnan');
    node_MSF(n)      = mean(Group_MSF(edge_from_n & valid_edges), 'omitnan');
    node_context(n)  = mean(C_edge(edge_from_n & valid_edges), 'omitnan');
end

% Correlations at node level
[r_node_MSF, p_node_MSF] = corr(node_MSF, node_duration, 'Type', 'Spearman');
[r_node_ctx, p_node_ctx] = corr(node_context, node_duration, 'Type', 'Spearman');

fprintf('Node-level correlations:\n');
fprintf('  Node MSF vs Duration:      rho = %+.3f, p = %.2e\n', r_node_MSF, p_node_MSF);
fprintf('  Node Context vs Duration:  rho = %+.3f, p = %.2e\n', r_node_ctx, p_node_ctx);

%% ========================================================================
%  SECTION 6: Publication-STYLE FIGURE
% ========================================================================

% --- Colors (NPG-inspired) ---
col_A    = [0.23, 0.33, 0.53];   % Navy - MSF
col_B    = [0.89, 0.29, 0.20];   % Coral - Context
col_C    = [0.00, 0.62, 0.45];   % Teal - Event count
col_gray = [0.55, 0.55, 0.55];

fig = figure('Units', 'normalized', 'Position', [0.05, 0.10, 0.92, 0.70], ...
             'Color', 'w', 'Name', 'Bridging_MSF_EventDuration');
set(gcf, 'DefaultAxesFontName', 'Helvetica', ...
         'DefaultTextFontName', 'Helvetica');

% -----------------------------------------------------------------------
%  PANEL A: MSF vs Event Duration (edge-level)
% -----------------------------------------------------------------------
axA = subplot(2, 3, 1); hold(axA, 'on');

% Density-based coloring
n_bins = 40;
x_edges = linspace(min(MSF_valid), max(MSF_valid), n_bins + 1);
y_edges = linspace(min(Dur_valid), max(Dur_valid), n_bins + 1);
[counts, ~, ~] = histcounts2(MSF_valid, Dur_valid, x_edges, y_edges);
[~, ~, bx] = histcounts(MSF_valid, x_edges);
[~, ~, by] = histcounts(Dur_valid, y_edges);
density = counts(sub2ind([n_bins, n_bins], bx, by));

[~, sort_idx] = sort(density, 'ascend');
scatter(axA, MSF_valid(sort_idx), Dur_valid(sort_idx), ...
        15, density(sort_idx), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
colormap(axA, turbo);
cbar = colorbar(axA);
cbar.Label.String = 'Density';
cbar.Label.FontSize = 8;
cbar.Location = 'eastoutside';

% Robust fit
mdl_A = fitlm(MSF_valid, Dur_valid, 'RobustOpts', 'on');
x_line = linspace(min(MSF_valid), max(MSF_valid), 100);
plot(axA, x_line, predict(mdl_A, x_line'), '-k', 'LineWidth', 2.5);

% Statistical annotation
txt_A = sprintf('\\rho = %+.3f\np < 0.001\nN = %d', rho_A, n_valid);
text(axA, min(MSF_valid) + 0.05*range(MSF_valid), ...
     max(Dur_valid) - 0.10*range(Dur_valid), txt_A, ...
     'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
     'Margin', 4);

xlabel(axA, 'MSF (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(axA, 'Mean event duration (s)', 'FontSize', 11, 'FontWeight', 'bold');
title(axA, 'A | Spectral vs. event-based', 'FontSize', 12, 'FontWeight', 'bold');
set(axA, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
grid(axA, 'off');

% -----------------------------------------------------------------------
%  PANEL B: Structural Context vs Event Duration
% -----------------------------------------------------------------------
axB = subplot(2, 3, 2); hold(axB, 'on');

x_edges = linspace(min(C_valid), max(C_valid), n_bins + 1);
y_edges = linspace(min(Dur_valid), max(Dur_valid), n_bins + 1);
[counts, ~, ~] = histcounts2(C_valid, Dur_valid, x_edges, y_edges);
[~, ~, bx] = histcounts(C_valid, x_edges);
[~, ~, by] = histcounts(Dur_valid, y_edges);
density = counts(sub2ind([n_bins, n_bins], bx, by));

[~, sort_idx] = sort(density, 'ascend');
scatter(axB, C_valid(sort_idx), Dur_valid(sort_idx), ...
        15, density(sort_idx), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
colormap(axB, turbo);
cbar = colorbar(axB);
cbar.Label.String = 'Density';
cbar.Label.FontSize = 8;
cbar.Location = 'eastoutside';

mdl_B = fitlm(C_valid, Dur_valid, 'RobustOpts', 'on');
x_line = linspace(min(C_valid), max(C_valid), 100);
plot(axB, x_line, predict(mdl_B, x_line'), '-k', 'LineWidth', 2.5);

txt_B = sprintf('\\rho = %+.3f\np < 0.001\nN = %d', rho_B, n_valid);
text(axB, max(C_valid) - 0.35*range(C_valid), ...
     max(Dur_valid) - 0.10*range(Dur_valid), txt_B, ...
     'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
     'Margin', 4);

xlabel(axB, 'Structural edge context', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(axB, 'Mean event duration (s)', 'FontSize', 11, 'FontWeight', 'bold');
title(axB, 'B | Structural context vs. event', 'FontSize', 12, 'FontWeight', 'bold');
set(axB, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
grid(axB, 'off');

% -----------------------------------------------------------------------
%  PANEL C: Event Count vs Event Duration (control)
% -----------------------------------------------------------------------
axC = subplot(2, 3, 3); hold(axC, 'on');

x_edges = linspace(min(Count_valid), max(Count_valid), n_bins + 1);
y_edges = linspace(min(Dur_valid), max(Dur_valid), n_bins + 1);
[counts, ~, ~] = histcounts2(Count_valid, Dur_valid, x_edges, y_edges);
[~, ~, bx] = histcounts(Count_valid, x_edges);
[~, ~, by] = histcounts(Dur_valid, y_edges);
density = counts(sub2ind([n_bins, n_bins], bx, by));

[~, sort_idx] = sort(density, 'ascend');
scatter(axC, Count_valid(sort_idx), Dur_valid(sort_idx), ...
        15, density(sort_idx), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
colormap(axC, turbo);
cbar = colorbar(axC);
cbar.Label.String = 'Density';
cbar.Label.FontSize = 8;
cbar.Location = 'eastoutside';

mdl_C = fitlm(Count_valid, Dur_valid, 'RobustOpts', 'on');
x_line = linspace(min(Count_valid), max(Count_valid), 100);
plot(axC, x_line, predict(mdl_C, x_line'), '-k', 'LineWidth', 2.5);

txt_C = sprintf('\\rho = %+.3f\np < 0.001\nN = %d', rho_C, n_valid);
text(axC, min(Count_valid) + 0.05*range(Count_valid), ...
     max(Dur_valid) - 0.10*range(Dur_valid), txt_C, ...
     'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
     'Margin', 4);

xlabel(axC, 'Event count', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(axC, 'Mean event duration (s)', 'FontSize', 11, 'FontWeight', 'bold');
title(axC, 'C | Event count (control)', 'FontSize', 12, 'FontWeight', 'bold');
set(axC, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
grid(axC, 'off');

% -----------------------------------------------------------------------
%  PANEL D: Bootstrap distributions
% -----------------------------------------------------------------------
axD = subplot(2, 3, 4); hold(axD, 'on');

histogram(axD, boot_A, 40, 'FaceColor', col_A, 'EdgeColor', 'none', ...
          'FaceAlpha', 0.6, 'Normalization', 'pdf');
histogram(axD, boot_B, 40, 'FaceColor', col_B, 'EdgeColor', 'none', ...
          'FaceAlpha', 0.6, 'Normalization', 'pdf');

xline(axD, rho_A, '-', 'Color', col_A, 'LineWidth', 2.5);
xline(axD, rho_B, '-', 'Color', col_B, 'LineWidth', 2.5);
xline(axD, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

xlabel(axD, 'Spearman \rho', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(axD, 'Probability density', 'FontSize', 11, 'FontWeight', 'bold');
title(axD, 'D | Bootstrap distributions', 'FontSize', 12, 'FontWeight', 'bold');

legend(axD, {sprintf('MSF vs Duration (\\rho = %+.2f)', rho_A), ...
             sprintf('Context vs Duration (\\rho = %+.2f)', rho_B)}, ...
             'Location', 'northwest', 'Box', 'off', 'FontSize', 9);

set(axD, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
grid(axD, 'off');

% -----------------------------------------------------------------------
%  PANEL E: Node-level comparison
% -----------------------------------------------------------------------
axE = subplot(2, 3, 5); hold(axE, 'on');

scatter(axE, node_MSF, node_duration, 60, col_A, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');

mdl_E = fitlm(node_MSF, node_duration, 'RobustOpts', 'on');
x_line = linspace(min(node_MSF), max(node_MSF), 50);
plot(axE, x_line, predict(mdl_E, x_line'), '-k', 'LineWidth', 2.5);

txt_E = sprintf('\\rho = %+.3f\np = %.2e', r_node_MSF, p_node_MSF);
text(axE, min(node_MSF) + 0.05*range(node_MSF), ...
     max(node_duration) - 0.10*range(node_duration), txt_E, ...
     'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
     'Margin', 4);

xlabel(axE, 'Node-level MSF', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(axE, 'Node-level duration (s)', 'FontSize', 11, 'FontWeight', 'bold');
title(axE, 'E | Node-level convergence', 'FontSize', 12, 'FontWeight', 'bold');
set(axE, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
grid(axE, 'off');

% -----------------------------------------------------------------------
%  PANEL F: Summary schematic bar plot
% -----------------------------------------------------------------------
axF = subplot(2, 3, 6); hold(axF, 'on');

% Summary metrics
metrics = [rho_A, rho_B, rho_partial, r_node_MSF, r_node_ctx];
labels  = {'MSF vs Dur', 'Ctx vs Dur', 'Partial\n(Ctx, Dur | MSF)', ...
           'Node MSF vs Dur', 'Node Ctx vs Dur'};
colors  = [col_A; col_B; col_C; col_A*0.7; col_B*0.7];

b = bar(axF, 1:5, metrics, 0.65, 'FaceColor', 'flat', ...
        'EdgeColor', 'k', 'LineWidth', 1.2);
for k = 1:5
    b.CData(k, :) = colors(k, :);
end

% Add numeric labels
for k = 1:5
    text(axF, k, metrics(k) + 0.03*sign(metrics(k)), ...
         sprintf('%.3f', metrics(k)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', 'bold');
end

yline(axF, 0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
set(axF, 'XTick', 1:5, 'XTickLabel', labels, ...
        'FontSize', 9, 'TickDir', 'out', 'Box', 'off');
xtickangle(axF, 15);
ylabel(axF, 'Spearman \rho', 'FontSize', 11, 'FontWeight', 'bold');
title(axF, 'F | Summary of convergence', 'FontSize', 12, 'FontWeight', 'bold');
ylim(axF, [-0.7, 0.1]);
grid(axF, 'off');

% --- Global title ---
sgtitle('Convergent evidence for edge-level timescale organization', ...
        'FontSize', 14, 'FontWeight', 'bold');

% --- Save ---
saveas(fig, 'Figure_SX_Bridging_MSF_EventDuration.png');
saveas(fig, 'Figure_SX_Bridging_MSF_EventDuration.svg');
try
    exportgraphics(fig, 'Figure_SX_Bridging_MSF_EventDuration.eps', ...
                   'ContentType', 'vector');
catch
    print(fig, 'Figure_SX_Bridging_MSF_EventDuration.eps', '-depsc2', '-r300');
end

fprintf('\nFigure saved: Figure_SX_Bridging_MSF_EventDuration.{png,svg,eps}\n');

%% ========================================================================
%  SECTION 7: SAVE NUMERICAL SUMMARY
% ========================================================================
T = table( ...
    {'MSF vs Duration'; 'Context vs Duration'; 'Partial (Ctx, Dur | MSF)'; ...
     'Direct edges only'; 'Non-direct edges only'; ...
     'Node-level MSF vs Duration'; 'Node-level Context vs Duration'}, ...
    [rho_A; rho_B; rho_partial; rho_dir; rho_non; r_node_MSF; r_node_ctx], ...
    [p_A; p_B; p_partial; p_dir; p_non; p_node_MSF; p_node_ctx], ...
    'VariableNames', {'Analysis', 'Spearman_rho', 'p_value'});

disp(T);
writetable(T, 'Bridging_MSF_EventDuration_Summary.csv');
fprintf('Summary saved: Bridging_MSF_EventDuration_Summary.csv\n');

%% ========================================================================
%  SECTION 8: FINAL SUMMARY
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  SUMMARY FOR MANUSCRIPT\n');
fprintf('===============================================================\n');
fprintf('\n(A) MSF vs Event Duration (edge-level):\n');
fprintf('    rho = %+.3f, p = %.2e, 95%% CI = [%+.3f, %+.3f]\n', ...
        rho_A, p_A, ci_A(1), ci_A(2));
fprintf('\n(B) Structural Context vs Event Duration (edge-level):\n');
fprintf('    rho = %+.3f, p = %.2e, 95%% CI = [%+.3f, %+.3f]\n', ...
        rho_B, p_B, ci_B(1), ci_B(2));
fprintf('\n(C) Partial correlation (Context vs Duration | MSF):\n');
fprintf('    rho_partial = %+.3f, p = %.2e\n', rho_partial, p_partial);
fprintf('\n(D) Node-level correlations:\n');
fprintf('    MSF vs Duration:     rho = %+.3f, p = %.2e\n', r_node_MSF, p_node_MSF);
fprintf('    Context vs Duration: rho = %+.3f, p = %.2e\n', r_node_ctx, p_node_ctx);
fprintf('\n===============================================================\n');
fprintf('  INTERPRETATION\n');
fprintf('===============================================================\n');
if rho_A < -0.4 && p_A < 0.001
    fprintf('✓ Strong convergence: MSF and event duration capture\n');
    fprintf('  complementary aspects of the same timescale property.\n');
end
if rho_B < -0.2 && p_B < 0.001
    fprintf('✓ Structural context predicts event duration,\n');
    fprintf('  linking spectral and event-based views.\n');
end
if abs(rho_partial) < abs(rho_B) * 0.7
    fprintf('! Partial correlation is attenuated, suggesting\n');
    fprintf('  some shared variance between MSF and duration.\n');
end
fprintf('===============================================================\n');
% --- Save duration metrics from the previous edge-event analysis ---
save('EdgeEvent_Duration_Results.mat', ...
     'Duration_group_mean', 'EventCount_group_mean', ...
     'Duration_per_subj', 'EventCount_per_subj', '-v7.3');
%% ========================================================================
%  DEGREE-PRESERVING NULL MODEL & STRUCTURE-TIMESCALE COUPLING ANALYSIS
%  Method: Maslov-Sneppen Edge Rewiring on Group SC
% ========================================================================

clearvars -except SC edge_median_freq; clc; close all;
fprintf('Loading structural connectivity...\n');
Results = load('EdgeSpace_Structure_Function_results_68regions.mat');
edge_median_freq = Results.Group_MSF; 
load('SC_SIFT2_68_100subj.mat');

% SC dimensions: N x N x S
SC = SC_all;
% --- Step 1: Data Preparation & Group-Average SC ---
fprintf('Step 1: Computing Group-Average SC Matrix...\n');

if ndims(SC) == 3
    Group_SC = mean(SC, 3); % Group-average across subjects (3rd dimension)
else
    Group_SC = SC;
end

N = size(Group_SC, 1);
Group_SC(logical(eye(N))) = 0; % Zero diagonal

% Extract lower triangular indices (unique edge pairs)
tril_mask = tril(true(N), -1);
[i_idx, j_idx] = find(tril_mask);
N_edges = length(i_idx);
tril_indices = sub2ind([N, N], i_idx, j_idx);

% Extract edge-wise vectors
Group_SC_vec = Group_SC(tril_indices);
mask_direct = (Group_SC_vec > 0);        % Direct structural edges (w_ij > 0)
mask_nondirect = (Group_SC_vec == 0);     % Non-direct structural pairs (w_ij = 0)

% --- Step 2: Empirical Structural Context Calculation ---
% C_ij^S = s_i + s_j - w_ij
s_emp = sum(Group_SC, 2);                     % Node strengths (N x 1)
S_mat_emp = s_emp + s_emp';                   % Matrix of s_i + s_j
Context_mat_emp = S_mat_emp - Group_SC;       % Empirical Context Matrix
C_emp = Context_mat_emp(tril_indices);        % Vector of empirical C_ij^S

% Compute Empirical Spearman Correlations
rho_emp_all = corr(C_emp, edge_median_freq, 'Type', 'Spearman', 'Rows', 'complete');
rho_emp_direct = corr(C_emp(mask_direct), edge_median_freq(mask_direct), 'Type', 'Spearman', 'Rows', 'complete');
rho_emp_nondirect = corr(C_emp(mask_nondirect), edge_median_freq(mask_nondirect), 'Type', 'Spearman', 'Rows', 'complete');

% --- Step 3: Null Model Iterations ---
n_nulls = 10000;
null_rhos_all = zeros(n_nulls, 1);
null_rhos_direct = zeros(n_nulls, 1);
null_rhos_nondirect = zeros(n_nulls, 1);

% Matrix to store null structural contexts for mean scatter plotting
null_C_matrix = zeros(N_edges, n_nulls);

fprintf('Step 2: Running %d Degree-Preserving Null Iterations...\n', n_nulls);
tic;
for i = 1:n_nulls
    if mod(i, 100) == 0 || i == n_nulls
        fprintf('  Processing null iteration %d / %d (%.1f sec elapsed)...\n', i, n_nulls, toc);
    end
    
    % 1. Rewire Group SC preserving degree distribution
    SC_null = degree_preserving_rewiring(Group_SC, 10);
    
    % 2. Fast Vectorized Structural Context Calculation
    s_null = sum(SC_null, 2);
    Context_mat_null = (s_null + s_null') - SC_null;
    C_null = Context_mat_null(tril_indices);
    null_C_matrix(:, i) = C_null;
    
    % 3. Calculate Null Correlations with Empirical Timescales (MSF)
    null_rhos_all(i) = corr(C_null, edge_median_freq, 'Type', 'Spearman', 'Rows', 'complete');
    null_rhos_direct(i) = corr(C_null(mask_direct), edge_median_freq(mask_direct), 'Type', 'Spearman', 'Rows', 'complete');
    null_rhos_nondirect(i) = corr(C_null(mask_nondirect), edge_median_freq(mask_nondirect), 'Type', 'Spearman', 'Rows', 'complete');
end

% Mean Null Context across iterations
C_null_mean = mean(null_C_matrix, 2);

% --- Step 4: Non-parametric p-values & Effect Sizes ---
p_null_all = sum(null_rhos_all <= rho_emp_all) / n_nulls;
p_null_direct = sum(null_rhos_direct <= rho_emp_direct) / n_nulls;
p_null_nondirect = sum(null_rhos_nondirect <= rho_emp_nondirect) / n_nulls;

z_score_all = (rho_emp_all - mean(null_rhos_all)) / std(null_rhos_all);
z_score_direct = (rho_emp_direct - mean(null_rhos_direct)) / std(null_rhos_direct);
z_score_nondirect = (rho_emp_nondirect - mean(null_rhos_nondirect)) / std(null_rhos_nondirect);

% Print Summary Results to Command Window
fprintf('\n================ ANALYSIS RESULTS SUMMARY ================\n');
fprintf('ALL EDGES (N = %d):\n', N_edges);
fprintf('  Empirical rho: %.4f | Null mean: %.4f (SD: %.4f) | Z: %.2f | p_null: %.4f\n', ...
    rho_emp_all, mean(null_rhos_all), std(null_rhos_all), z_score_all, p_null_all);

fprintf('DIRECT EDGES (w_ij > 0, N = %d):\n', sum(mask_direct));
fprintf('  Empirical rho: %.4f | Null mean: %.4f (SD: %.4f) | Z: %.2f | p_null: %.4f\n', ...
    rho_emp_direct, mean(null_rhos_direct), std(null_rhos_direct), z_score_direct, p_null_direct);

fprintf('NON-DIRECT PAIRS (w_ij = 0, N = %d):\n', sum(mask_nondirect));
fprintf('  Empirical rho: %.4f | Null mean: %.4f (SD: %.4f) | Z: %.2f | p_null: %.4f\n', ...
    rho_emp_nondirect, mean(null_rhos_nondirect), std(null_rhos_nondirect), z_score_nondirect, p_null_nondirect);
fprintf('==========================================================\n\n');


%% ========================================================================
%  VISUALIZATION: Publication-Ready 4-Panel Nature-Style Figure
% ========================================================================

fig = figure('Position', [100, 100, 1250, 850], 'Color', 'w', 'Name', 'Degree-Preserving Null Model');
set(gcf, 'DefaultAxesFontName', 'Helvetica', 'DefaultTextFontName', 'Helvetica');

% Palette definition
color_emp_all = [0.15, 0.35, 0.55];      % Deep Navy/Slate
color_null = [0.55, 0.45, 0.65];         % Muted Purple
color_direct = [0.85, 0.32, 0.10];       % Coral/Red
color_nondirect = [0.12, 0.47, 0.71];    % Blue

% -------------------------------------------------------------------------
% PANEL A: All Edges Null Distribution
% -------------------------------------------------------------------------
subplot(2,2,1);
h1 = histogram(null_rhos_all, 'Normalization', 'pdf', 'FaceColor', color_null, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
[f_a, x_a] = ksdensity(null_rhos_all);
plot(x_a, f_a, 'k-', 'LineWidth', 2);
xline(rho_emp_all, 'r-', 'LineWidth', 2.5);

xlabel('Spearman \rho (Structural Context vs. MSF)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 10, 'FontWeight', 'bold');
title('A. Overall Coupling: Null Model Test (All Pairs)', 'FontSize', 12, 'FontWeight', 'bold');

txt_A = sprintf('All Edges (N = %d)\nEmpirical \\rho = %.3f\nNull mean = %.3f\nEffect Size (Z) = %.2f\np_{null} %s', ...
    N_edges, rho_emp_all, mean(null_rhos_all), z_score_all, format_p_val(p_null_all));
text(min(xlim)+0.04*range(xlim), max(ylim)*0.82, txt_A, 'FontSize', 9.5, ...
    'EdgeColor', [0.8 0.8 0.8], 'BackgroundColor', [0.98 0.98 0.98], 'Margin', 6);
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 10);

% -------------------------------------------------------------------------
% PANEL B: Empirical vs. Null Structural Context Scatter
% -------------------------------------------------------------------------
subplot(2,2,2);
% Downsample points for visual clarity if N_edges > 5000
step_ds = max(1, round(N_edges / 3000));
idx_plot = 1:step_ds:N_edges;

scatter(C_emp(idx_plot), edge_median_freq(idx_plot), 15, color_emp_all, ...
    'filled', 'MarkerFaceAlpha', 0.25); hold on;
scatter(C_null_mean(idx_plot), edge_median_freq(idx_plot), 15, [0.6 0.6 0.6], ...
    'filled', 'MarkerFaceAlpha', 0.20);

% Polyfit regression trends
p_emp_all = polyfit(C_emp, edge_median_freq, 1);
p_null_all_fit = polyfit(C_null_mean, edge_median_freq, 1);

x_eval = [min(C_emp), max(C_emp)];
plot(x_eval, polyval(p_emp_all, x_eval), 'Color', color_emp_all, 'LineWidth', 2.5);
plot(x_eval, polyval(p_null_all_fit, x_eval), 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'LineStyle', '--');

xlabel('Structural Context (C_{ij}^S)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Functional Timescale (MSF)', 'FontSize', 10, 'FontWeight', 'bold');
title('B. Empirical vs. Rewired Topology Dynamics', 'FontSize', 12, 'FontWeight', 'bold');
legend({sprintf('Empirical (\\rho = %.2f)', rho_emp_all), ...
        sprintf('Null Average (\\rho = %.2f)', mean(null_rhos_all))}, ...
        'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 10);

% -------------------------------------------------------------------------
% PANEL C: Direct Structural Edges (w_ij > 0) Null Test
% -------------------------------------------------------------------------
subplot(2,2,3);
histogram(null_rhos_direct, 'Normalization', 'pdf', 'FaceColor', color_direct, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
[f_c, x_c] = ksdensity(null_rhos_direct);
plot(x_c, f_c, 'k-', 'LineWidth', 2);
xline(rho_emp_direct, 'r-', 'LineWidth', 2.5);

xlabel('Spearman \rho (Structural Context vs. MSF)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 10, 'FontWeight', 'bold');
title('C. Direct Structural Edges (w_{ij} > 0)', 'FontSize', 12, 'FontWeight', 'bold');

txt_C = sprintf('Direct Edges (N = %d)\nEmpirical \\rho = %.3f\nNull mean = %.3f\nEffect Size (Z) = %.2f\np_{null} %s', ...
    sum(mask_direct), rho_emp_direct, mean(null_rhos_direct), z_score_direct, format_p_val(p_null_direct));
text(min(xlim)+0.04*range(xlim), max(ylim)*0.82, txt_C, 'FontSize', 9.5, ...
    'EdgeColor', [0.8 0.8 0.8], 'BackgroundColor', [0.98 0.98 0.98], 'Margin', 6);
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 10);

% -------------------------------------------------------------------------
% PANEL D: Non-Direct Structural Pairs (w_ij = 0) Null Test
% -------------------------------------------------------------------------
subplot(2,2,4);
histogram(null_rhos_nondirect, 'Normalization', 'pdf', 'FaceColor', color_nondirect, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
[f_d, x_d] = ksdensity(null_rhos_nondirect);
plot(x_d, f_d, 'k-', 'LineWidth', 2);
xline(rho_emp_nondirect, 'r-', 'LineWidth', 2.5);

xlabel('Spearman \rho (Structural Context vs. MSF)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 10, 'FontWeight', 'bold');
title('D. Non-Direct Structural Pairs (w_{ij} = 0)', 'FontSize', 12, 'FontWeight', 'bold');

txt_D = sprintf('Non-Direct Pairs (N = %d)\nEmpirical \\rho = %.3f\nNull mean = %.3f\nEffect Size (Z) = %.2f\np_{null} %s', ...
    sum(mask_nondirect), rho_emp_nondirect, mean(null_rhos_nondirect), z_score_nondirect, format_p_val(p_null_nondirect));
text(min(xlim)+0.04*range(xlim), max(ylim)*0.82, txt_D, 'FontSize', 9.5, ...
    'EdgeColor', [0.8 0.8 0.8], 'BackgroundColor', [0.98 0.98 0.98], 'Margin', 6);
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 10);

sgtitle('Structural Context vs. Functional Timescale: Degree-Preserving Null Verification', ...
    'FontSize', 14, 'FontWeight', 'bold');


%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function W_null = degree_preserving_rewiring(W, ITER)
    % Maslov-Sneppen rewiring algorithm preserving node degree distribution
    W_null = W;
    [i_edges, j_edges] = find(tril(W_null, -1));
    K = length(i_edges);
    max_iter = round(ITER * K);
    
    for iter = 1:max_iter
        e1 = randi(K);
        e2 = randi(K);
        
        a = i_edges(e1); b = j_edges(e1);
        c = i_edges(e2); d = j_edges(e2);
        
        if length(unique([a b c d])) == 4
            if W_null(a,c) == 0 && W_null(b,d) == 0
                weight_ab = W_null(a,b);
                weight_cd = W_null(c,d);
                
                W_null(a,c) = weight_ab; W_null(c,a) = weight_ab;
                W_null(b,d) = weight_cd; W_null(d,b) = weight_cd;
                
                W_null(a,b) = 0; W_null(b,a) = 0;
                W_null(c,d) = 0; W_null(d,c) = 0;
                
                i_edges(e1) = max(a,c); j_edges(e1) = min(a,c);
                i_edges(e2) = max(b,d); j_edges(e2) = min(b,d);
            end
        end
    end
end

function str = format_p_val(p)
    if p == 0
        str = '< 0.0001';
    else
        str = sprintf('= %.4f', p);
    end
end
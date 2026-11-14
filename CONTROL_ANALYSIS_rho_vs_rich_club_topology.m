%% ========================================================================
%  RICH-CLUB ANALYSIS: Does hub organization explain edge-level coupling?
%  Two analyses:
%    (A) Subject-level: ρ_s vs rich-club coefficient
%    (B) Edge-level: partition edges into Rich-club / Feeder / Local
% ========================================================================

clearvars; clc; close all;

%% ========================================================================
%  SECTION 1: LOAD DATA
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  RICH-CLUB ANALYSIS OF EDGE-LEVEL COUPLING\n');
fprintf('===============================================================\n');

load('SC_SIFT2_68_100subj.mat');   % SC_all: Nnode x Nnode x 100
D = load('EdgeSpace_Structure_Function_results.mat');
C_edge      = D.Structural_Edge_Context(:);   
Subject_MSF = D.Subject_MSF;                  
mask_direct = logical(D.Direct_mask(:));
mask_nondir = logical(D.NonDirect_mask(:));

[N_subj, N_edges] = size(Subject_MSF);
N_nodes = size(SC_all, 1);

% Rebuild edge indices (lower triangle)
tril_mask = tril(true(N_nodes), -1);
[i_idx, j_idx] = find(tril_mask);

% Subject-level ρ_s
rho_s = zeros(N_subj, 1);
for s = 1:N_subj
    rho_s(s) = corr(C_edge, Subject_MSF(s, :)', 'Type', 'Spearman');
end
fprintf('ρ_s: mean = %.4f, SD = %.4f\n', mean(rho_s), std(rho_s));

%% ========================================================================
%  SECTION 2: COMPUTE RICH-CLUB COEFFICIENT PER SUBJECT
% ========================================================================
fprintf('\n--- Computing rich-club metrics per subject ---\n');

% Define k range for rich-club analysis (top-k nodes by degree)
% k_range = 10:5:60;   % From 10 to 60 hubs (Glasser has 360 nodes)
k_range = 5:5:40;
n_k = length(k_range);

% Pre-allocate
Phi_norm_all        = zeros(N_subj, 1);   % Normalized rich-club coefficient
Phi_raw_all         = zeros(N_subj, n_k); % Raw Φ(k)
rich_club_strength  = zeros(N_subj, 1);   % Mean of Φ_norm across k

for s = 1:N_subj
    A = SC_all(:, :, s);
    A = A - diag(diag(A));
    A(A < 0) = 0;
    
    % Binarize for topological rich-club analysis
    % (keep edges above a threshold to remove weak connections)
    A_bin = double(A > 0);
    
    % Compute rich-club coefficient Φ(k)
    Phi_raw = zeros(1, n_k);
    for ki = 1:n_k
        k = k_range(ki);
        
        % Rank nodes by degree
        deg = sum(A_bin, 2);
        [~, sorted_nodes] = sort(deg, 'descend');
        hub_nodes = sorted_nodes(1:k);
        
        % Subgraph of top-k nodes
        A_sub = A_bin(hub_nodes, hub_nodes);
        E_hub = sum(A_sub(:)) / 2;   % Number of edges in subgraph
        E_max = k * (k - 1) / 2;     % Max possible edges
        
        Phi_raw(ki) = E_hub / E_max;
    end
    
    % Random null model: compare to randomized networks with same degree
    Phi_rand = zeros(1, n_k);
    n_rand = 20;   % Number of randomizations
    for r = 1:n_rand
        A_rand = randomize_edges(A_bin, 5 * sum(A_bin(:)));
        deg_r = sum(A_rand, 2);
        [~, sorted_r] = sort(deg_r, 'descend');
        
        for ki = 1:n_k
            k = k_range(ki);
            hub_r = sorted_r(1:k);
            A_sub_r = A_rand(hub_r, hub_r);
            Phi_rand(ki) = Phi_rand(ki) + ...
                (sum(A_sub_r(:))/2) / (k*(k-1)/2);
        end
    end
    Phi_rand = Phi_rand / n_rand;
    
    % Normalized rich-club coefficient
    Phi_norm = Phi_raw ./ (Phi_rand + eps);
    Phi_raw_all(s, :) = Phi_raw;
    
    % Aggregate metric: mean of normalized Φ across k range
    % (focus on k where Φ_norm > 1, indicating true rich-club)
    valid_k = Phi_norm > 1;
    if sum(valid_k) > 0
        rich_club_strength(s) = mean(Phi_norm(valid_k));
    else
        rich_club_strength(s) = 1;
    end
    
    Phi_norm_all(s) = mean(Phi_norm);
end

fprintf('Rich-club metrics computed for %d subjects.\n', N_subj);
fprintf('  Mean Φ_norm = %.3f ± %.3f\n', mean(Phi_norm_all), std(Phi_norm_all));

%% ========================================================================
%  SECTION 3: SUBJECT-LEVEL CORRELATION ρ_s vs RICH-CLUB
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  SUBJECT-LEVEL: ρ_s vs RICH-CLUB\n');
fprintf('===============================================================\n');

[r_rc, p_rc] = corr(rho_s, rich_club_strength, 'Type', 'Spearman');
[r_rc_mean, p_rc_mean] = corr(rho_s, Phi_norm_all, 'Type', 'Spearman');

fprintf('ρ_s vs Rich-club strength:  ρ = %+.3f, p = %.4f, R² = %.3f\n', ...
        r_rc, p_rc, r_rc^2);
fprintf('ρ_s vs Mean Φ_norm:         ρ = %+.3f, p = %.4f, R² = %.3f\n', ...
        r_rc_mean, p_rc_mean, r_rc_mean^2);

% Interpretation
if p_rc < 0.05
    fprintf('\n*** POSITIVE FINDING: ρ_s is significantly associated with rich-club.\n');
    fprintf('*** Interpretation: Subjects with stronger hub organization exhibit\n');
    fprintf('*** stronger edge-level coupling. This is a NEW DISCOVERY.\n');
else
    fprintf('\n*** NULL FINDING: ρ_s is independent of rich-club topology.\n');
    fprintf('*** Interpretation: Edge-level coupling is distributed across the\n');
    fprintf('*** network and not driven by hub topology.\n');
end

%% ========================================================================
%  SECTION 4: EDGE-LEVEL ANALYSIS - Rich/Feeder/Local Classification
% ========================================================================
fprintf('\n===============================================================\n');
fprintf('  EDGE-LEVEL: Rich-club vs Feeder vs Local connections\n');
fprintf('===============================================================\n');

% Use group-average connectome for hub classification
A_group = mean(SC_all, 3);
A_group = A_group - diag(diag(A_group));
A_group(A_group < 0) = 0;

% Define hubs as top 15% of nodes by degree
deg_group = sum(A_group, 2);
hub_threshold = prctile(deg_group, 85);   % top 15%
is_hub = deg_group >= hub_threshold;

fprintf('Hub threshold: degree >= %.2f\n', hub_threshold);
fprintf('Number of hubs: %d/%d (%.1f%%)\n', ...
        sum(is_hub), N_nodes, 100*sum(is_hub)/N_nodes);

% Classify edges
edge_is_hub_i = is_hub(i_idx);
edge_is_hub_j = is_hub(j_idx);

rich_club_mask  = edge_is_hub_i & edge_is_hub_j;       % hub-hub
feeder_mask     = xor(edge_is_hub_i, edge_is_hub_j);   % hub-periphery
local_mask      = ~edge_is_hub_i & ~edge_is_hub_j;     % periphery-periphery

fprintf('\nEdge classification:\n');
fprintf('  Rich-club (hub-hub):      %d edges (%.1f%%)\n', ...
        sum(rich_club_mask), 100*sum(rich_club_mask)/N_edges);
fprintf('  Feeder (hub-periphery):   %d edges (%.1f%%)\n', ...
        sum(feeder_mask), 100*sum(feeder_mask)/N_edges);
fprintf('  Local (periphery-periph): %d edges (%.1f%%)\n', ...
        sum(local_mask), 100*sum(local_mask)/N_edges);

% Compute correlations for each edge class
C_rc = C_edge(rich_club_mask);
C_fd = C_edge(feeder_mask);
C_lc = C_edge(local_mask);

% Use group MSF
Group_MSF = mean(Subject_MSF, 1)';
MSF_rc = Group_MSF(rich_club_mask);
MSF_fd = Group_MSF(feeder_mask);
MSF_lc = Group_MSF(local_mask);

[rho_rc, p_rc_edges] = corr(C_rc, MSF_rc, 'Type', 'Spearman');
[rho_fd, p_fd_edges] = corr(C_fd, MSF_fd, 'Type', 'Spearman');
[rho_lc, p_lc_edges] = corr(C_lc, MSF_lc, 'Type', 'Spearman');

fprintf('\nCorrelations (Context vs MSF) by edge class:\n');
fprintf('  Rich-club:  ρ = %+.3f, p = %.3e, N = %d\n', rho_rc, p_rc_edges, sum(rich_club_mask));
fprintf('  Feeder:     ρ = %+.3f, p = %.3e, N = %d\n', rho_fd, p_fd_edges, sum(feeder_mask));
fprintf('  Local:      ρ = %+.3f, p = %.3e, N = %d\n', rho_lc, p_lc_edges, sum(local_mask));

% Test whether rich-club coupling is stronger than local coupling
% (using bootstrap)
fprintf('\nBootstrap test: Rich-club coupling vs Local coupling\n');
N_boot = 5000;
boot_diff = zeros(N_boot, 1);
rng(42);
for b = 1:N_boot
    idx_rc = randsample(length(C_rc), length(C_rc), true);
    idx_lc = randsample(length(C_lc), length(C_lc), true);
    r1 = corr(C_rc(idx_rc), MSF_rc(idx_rc), 'Type', 'Spearman');
    r2 = corr(C_lc(idx_lc), MSF_lc(idx_lc), 'Type', 'Spearman');
    boot_diff(b) = r1 - r2;
end
p_diff = 2 * min(mean(boot_diff <= 0), mean(boot_diff >= 0));
fprintf('  Δρ (rich - local) = %+.3f, p = %.4f\n', rho_rc - rho_lc, p_diff);

%% ========================================================================
%  SECTION 5: VISUALIZATION (Nature-style)
% ========================================================================

fig = figure('Units', 'normalized', 'Position', [0.05, 0.15, 0.90, 0.45], ...
             'Color', 'w');
set(gcf, 'DefaultAxesFontName', 'Helvetica', ...
         'DefaultTextFontName', 'Helvetica');

% --- Panel A: ρ_s vs Rich-club coefficient ---
axA = subplot(1, 3, 1); hold(axA, 'on');
scatter(axA, rich_club_strength, rho_s, 60, [0.23, 0.33, 0.53], ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');

% Robust fit
mdl_a = fitlm(rich_club_strength, rho_s, 'RobustOpts', 'on');
x_l = linspace(min(rich_club_strength), max(rich_club_strength), 50);
plot(axA, x_l, predict(mdl_a, x_l'), '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5);

text(axA, min(rich_club_strength) + 0.05*range(rich_club_strength), ...
     min(rho_s) + 0.1*range(rho_s), ...
     sprintf('\\rho = %+.3f\np = %.3f', r_rc, p_rc), ...
     'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], ...
     'Margin', 4);

xlabel(axA, 'Rich-club strength (\Phi_{norm})', 'FontSize', 11, 'FontWeight','bold');
ylabel(axA, '\rho_s', 'FontSize', 11, 'FontWeight','bold');
title(axA, 'A | \rho_s vs rich-club', 'FontSize', 12, 'FontWeight','bold');
set(axA, 'TickDir','out', 'Box','off', 'FontSize', 10);
grid(axA, 'off');

% --- Panel B: Edge classification scatter ---
axB = subplot(1, 3, 2); hold(axB, 'on');

scatter(axB, C_rc, MSF_rc, 10, [0.89, 0.29, 0.20], 'filled', ...
        'MarkerFaceAlpha', 0.4);
scatter(axB, C_fd, MSF_fd, 10, [0.85, 0.65, 0.15], 'filled', ...
        'MarkerFaceAlpha', 0.3);
scatter(axB, C_lc, MSF_lc, 10, [0.23, 0.33, 0.53], 'filled', ...
        'MarkerFaceAlpha', 0.3);

xlabel(axB, 'Structural edge context', 'FontSize', 11, 'FontWeight','bold');
ylabel(axB, 'MSF', 'FontSize', 11, 'FontWeight','bold');
title(axB, 'B | Rich / Feeder / Local', 'FontSize', 12, 'FontWeight','bold');
legend(axB, {sprintf('Rich (\\rho=%.2f)', rho_rc), ...
             sprintf('Feeder (\\rho=%.2f)', rho_fd), ...
             sprintf('Local (\\rho=%.2f)', rho_lc)}, ...
             'Location','northeast', 'Box','off', 'FontSize', 9);
set(axB, 'TickDir','out', 'Box','off', 'FontSize', 10);
grid(axB, 'off');

% --- Panel C: Bar comparison ---
axC = subplot(1, 3, 3); hold(axC, 'on');
rho_vals = [rho_rc, rho_fd, rho_lc];
b = bar(axC, 1:3, rho_vals, 0.65, 'FaceColor', 'flat', ...
        'EdgeColor', 'k', 'LineWidth', 1.2);
b.CData(1, :) = [0.89, 0.29, 0.20];
b.CData(2, :) = [0.85, 0.65, 0.15];
b.CData(3, :) = [0.23, 0.33, 0.53];

for i = 1:3
    text(axC, i, rho_vals(i) + 0.02, sprintf('%.3f', rho_vals(i)), ...
         'HorizontalAlignment','center', 'FontSize', 10, 'FontWeight','bold');
end

set(axC, 'XTick', 1:3, 'XTickLabel', {'Rich-club','Feeder','Local'}, ...
        'FontSize', 10, 'TickDir','out', 'Box','off');
ylabel(axC, 'Spearman \rho (Context vs. MSF)', 'FontSize', 11, 'FontWeight','bold');
title(axC, 'C | Coupling by edge class', 'FontSize', 12, 'FontWeight','bold');
yline(axC, 0, '-', 'Color', [0.4 0.4 0.4]);
grid(axC, 'off');

sgtitle('Rich-club analysis: is edge-level coupling hub-centric?', ...
        'FontSize', 13, 'FontWeight','bold');

% saveas(fig, 'Figure_RichClub_Analysis.png');
% saveas(fig, 'Figure_RichClub_Analysis.svg');
fprintf('\nFigure saved: Figure_RichClub_Analysis.png\n');

%% ========================================================================
%  LOCAL FUNCTION: RANDOMIZE EDGES (Maslov-Sneppen rewiring)
% ========================================================================
function A_rand = randomize_edges(A, n_swaps)
    % Degree-preserving edge-swapping (Maslov-Sneppen algorithm)
    A_rand = A;
    [i_edges, j_edges] = find(triu(A_rand));
    E = length(i_edges);
    
    for s = 1:n_swaps
        e1 = randi(E);
        e2 = randi(E);
        if e1 == e2, continue; end
        
        a = i_edges(e1); b = j_edges(e1);
        c = i_edges(e2); d = j_edges(e2);
        
        % Avoid self-loops and existing edges
        if a == c || a == d || b == c || b == d, continue; end
        if A_rand(a, d) || A_rand(b, c), continue; end
        
        % Swap
        A_rand(a, b) = 0; A_rand(b, a) = 0;
        A_rand(c, d) = 0; A_rand(d, c) = 0;
        A_rand(a, d) = 1; A_rand(d, a) = 1;
        A_rand(b, c) = 1; A_rand(c, b) = 1;
        
        % Update edge list
        i_edges(e1) = min(a,d); j_edges(e1) = max(a,d);
        i_edges(e2) = min(b,c); j_edges(e2) = max(b,c);
    end
end
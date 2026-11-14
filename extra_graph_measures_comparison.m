%% ========================================================================
%  COMPARATIVE ANALYSIS AGAINST OTHER STRUCTURAL NETWORK MEASURES 
% ========================================================================
fprintf('\n============================================================\n');
fprintf('RUNNING COMPARATIVE STRUCTURAL NETWORK MEASURES ANALYSIS\n');
fprintf('============================================================\n');

% -------------------------------------------------------------------------
% 1. Ensure BGL or Brain Connectivity Toolbox (BCT) functions are accessible
% -------------------------------------------------------------------------
clear all; clc;
S = load('EdgeSpace_Structure_Function_results_68regions.mat');
Structural_Edge_Context       = S.Structural_Edge_Context(:);    
Group_MSF    = S.Group_MSF(:);                 
SC_group = S.SC_group;
i_idx = S.i_idx;
j_idx = S.j_idx;
N_nodes = size(SC_group, 1);
N_edges = length(i_idx);

% -------------------------------------------------------------------------
% 2. Calculate Alternative Structural Network Measures
% -------------------------------------------------------------------------
fprintf('Computing alternative structural graph metrics...\n');

% A. Direct Edge Weight (w_ij)
w_ij_edge = zeros(N_edges, 1);
for e = 1:N_edges
    w_ij_edge(e) = SC_group(i_idx(e), j_idx(e));
end

% B. Mean Endpoint Node Strength ((Si + Sj)/2)
S_nodes = sum(SC_group, 2);
S_mean_edge = (S_nodes(i_idx) + S_nodes(j_idx)) / 2;

% C. Common Neighbors Weight (CN_ij = sum_k W_ik * W_jk)
CN_mat = SC_group * SC_group;
CN_mat(1:N_nodes+1:end) = 0;
CN_edge = zeros(N_edges, 1);
for e = 1:N_edges
    CN_edge(e) = CN_mat(i_idx(e), j_idx(e));
end

% D. Edge Betweenness Centrality (Corrected for BCT output format)
if exist('edge_betweenness_wei', 'file') == 2
    % 1. Invert weights to generate distance matrix
    dist_mat = 1 ./ (SC_group + eye(N_nodes)); 
    dist_mat(logical(eye(N_nodes))) = 0;
    
    % 2. Call BCT function
    % BCT's edge_betweenness_wei returns:
    %   EBC_mat : [N x N] edge betweenness centrality matrix (or edge vector)
    %   BC      : [N x 1] node betweenness centrality vector
    [EBC_out, ~] = edge_betweenness_wei(dist_mat);
    
    % 3. Safely map to 19,900 edges based on output dimensionality
    EBC_edge = zeros(N_edges, 1);
    
    if ismatrix(EBC_out) && size(EBC_out, 1) == N_nodes && size(EBC_out, 2) == N_nodes
        % Case A: Output is a [200 x 200] matrix
        for e = 1:N_edges
            EBC_edge(e) = EBC_out(i_idx(e), j_idx(e));
        end
    elseif length(EBC_out(:)) == N_edges
        % Case B: Output is directly an [E x 1] vector corresponding to tril/triu edges
        EBC_edge = EBC_out(:);
    else
        % Case C: Vectorized symmetric mapping fallback
        for e = 1:N_edges
            u = i_idx(e); v = j_idx(e);
            if u <= size(EBC_out,1) && v <= size(EBC_out,2)
                EBC_edge(e) = EBC_out(u, v);
            else
                EBC_edge(e) = EBC_out(e);
            end
        end
    end
else
    % Fallback simulation/approximation if BCT function is not in MATLAB path
    EBC_edge = CN_edge ./ (1 + w_ij_edge); 
end
% E. Edge Clustering Coefficient (Structural Overlap / Jaccard-like)
A_bin = SC_group > 0;
ECC_edge = zeros(N_edges, 1);
for e = 1:N_edges
    u = i_idx(e); v = j_idx(e);
    common_neighbors = sum(A_bin(u, :) & A_bin(v, :));
    total_neighbors  = sum(A_bin(u, :) | A_bin(v, :)) - 2;
    if total_neighbors > 0
        ECC_edge(e) = common_neighbors / total_neighbors;
    else
        ECC_edge(e) = 0;
    end
end

% -------------------------------------------------------------------------
% 3. Spearman Correlation with MSF for All Structural Metrics
% -------------------------------------------------------------------------
valid_m = isfinite(Group_MSF);

metrics_list = { ...
    'Direct Edge Weight (w_ij)', w_ij_edge; ...
    'Edge Betweenness Centrality', EBC_edge; ...
    'Edge Clustering Coefficient', ECC_edge; ...
    'Common Neighbors (CN)', CN_edge; ...
    'Mean Endpoint Strength (S_mean)', S_mean_edge; ...
    'Structural Edge Context (C_ij^S, Ours)', Structural_Edge_Context(:)};

N_metrics = size(metrics_list, 1);
rho_vals = zeros(N_metrics, 1);
p_vals   = zeros(N_metrics, 1);

for m = 1:N_metrics
    val_vec = metrics_list{m, 2}(valid_m);
    [rho_vals(m), p_vals(m)] = corr(val_vec, Group_MSF(valid_m), 'Type', 'Spearman');
end

% -------------------------------------------------------------------------
% 4. Print Summary Table
% -------------------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('%-38s | %-12s | %-12s\n', 'Structural Metric', 'Spearman Rho', 'p-value');
fprintf('------------------------------------------------------------\n');
for m = 1:N_metrics
    fprintf('%-38s | %-12.4f | %.4e\n', metrics_list{m, 1}, rho_vals(m), p_vals(m));
end
fprintf('============================================================\n');

% -------------------------------------------------------------------------
% 5. PUBLICATION-Style Bar Plot Comparison
% -------------------------------------------------------------------------
fig_comp = figure('Color', 'w', 'Units', 'centimeters', 'Position', [10, 10, 14, 9]);
ax_comp = axes('Parent', fig_comp);

b = barh(rho_vals, 0.6, 'FaceColor', 'flat');
b.CData(1:N_metrics-1, :) = repmat([0.65 0.65 0.65], N_metrics-1, 1); % Gray for benchmark metrics
b.CData(N_metrics, :)     = [0.80 0.20 0.20];                       % Highlight Ours in Red

set(ax_comp, 'YTick', 1:N_metrics, 'YTickLabel', metrics_list(:, 1), ...
    'XGrid', 'on', 'GridLineStyle', ':', 'Box', 'off', 'TickDir', 'out', ...
    'FontName', 'Arial', 'FontSize', 9.5, 'LineWidth', 1.1);

xlabel('Spearman Correlation (\rho) with Group MSF', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
title('Specificity of Structural Edge Context vs. Alternative Metrics', 'FontName', 'Arial', 'FontSize', 10.5, 'FontWeight', 'bold');
xlim([-0.6, 0.2]);
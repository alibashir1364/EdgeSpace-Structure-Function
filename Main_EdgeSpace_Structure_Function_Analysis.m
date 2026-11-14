%% ========================================================================
%  MAIN EDGE-SPACE STRUCTURE-FUNCTION ANALYSIS
%
%  Purpose:
%  Construct the structural edge-space graph and characterize the
%  temporal dynamics of all possible functional interactions.
%
%  Main correspondence:
%
%       Structural edge context C^S
%                    <-->
%                 Functional MSF
%
%  Dataset:
%       HCP 100 unrelated subjects
%
%  ------------------------------------------------------------------------
%  Structural pipeline:
%
%       Subject SC
%           ->
%       Group SC
%           ->
%       Complete edge space E*
%           ->
%       Edge-space graph G_E^S
%           ->
%       Structural edge context C^S
%
%  Functional pipeline:
%
%       BOLD
%           ->
%       z-scored BOLD
%           ->
%       ETS
%           ->
%       Welch PSD
%           ->
%       Subject MSF
%           ->
%       Group-average PSD
%           ->
%       Group MSF
%
%  ------------------------------------------------------------------------
%  The script stores variables required for:
%       1. Primary analysis
%       2. Network-aware statistical testing
%       3. Null models
%       4. Subject-level analysis
%       5. Sensitivity analysis
%
% ========================================================================

clear;
clc;
close all;

fprintf('\n');
fprintf('============================================================\n');
fprintf('EDGE-SPACE STRUCTURE-FUNCTION ANALYSIS\n');
fprintf('============================================================\n');


%% ========================================================================
% STEP 0 — USER SETTINGS
% ========================================================================

% -------------------------------------------------------------------------
% Dataset
% -------------------------------------------------------------------------

% fMRI_file = 'ROI_TS_shefer_200_HCP_100subj.mat';
fMRI_file = 'ROI_TS_68_HCP_100subj';
% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------

output_file = 'EdgeSpace_Structure_Function_results_68regions.mat';

% -------------------------------------------------------------------------
% Dataset parameters
% -------------------------------------------------------------------------

% N_nodes = 200;
N_nodes = 68;
N_subs  = 100;
N_time  = 1200;

TR = 0.72;
fs = 1 / TR;

% -------------------------------------------------------------------------
% Welch parameters
% -------------------------------------------------------------------------

window_length = 128;
window   = hanning(window_length);

noverlap = 64;
nfft     = 512;

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------

% Construct and store the explicit edge-space adjacency matrices.
%
% A_E_top and A_E_weighted are sparse matrices.
%
% For N = 200:
%
%       N_edges = 19,900
%       degree  = 396
%
store_edge_space_graph = true;

% Store subject-level MSF.
store_subject_level = true;

% Store subject-level total spectral energy.
store_subject_energy = true;

% Compute a sensitivity MSF after removing the DC component.
compute_noDC_sensitivity = true;


%% ========================================================================
% STEP 1 — LOAD DATA
% ========================================================================

fprintf('\nLoading fMRI data...\n');

load(fMRI_file);

fMRI = roiTS_all;

fprintf('Loading structural connectivity...\n');

%% ========================================================================
% STEP 2 — VERIFY DATA DIMENSIONS
% ========================================================================

[Nf, Tf, Sf] = size(fMRI);
Ns1 = Nf;
Ns2 = Nf;
Ss = Sf;

fprintf('\nData dimensions:\n');
fprintf('fMRI = [%d x %d x %d]\n', Nf, Tf, Sf);
fprintf('SC   = [%d x %d x %d]\n', Ns1, Ns2, Ss);

if Nf ~= N_nodes || Tf ~= N_time || Sf ~= N_subs
    error('Unexpected fMRI dimensions.');
end

if Ns1 ~= N_nodes || Ns2 ~= N_nodes || Ss ~= N_subs
    error('Unexpected SC dimensions.');
end

fprintf('Data dimensions verified successfully.\n');


%% ========================================================================
% STEP 3 — COMPLETE EDGE SPACE
%
% E* = { e_ij : 1 <= i < j <= N }
%
% Number of possible undirected edges:
%
%       N_edges = N(N-1)/2
%
% ========================================================================

fprintf('\nConstructing complete edge space...\n');

[i_idx, j_idx] = find(tril(ones(N_nodes), -1));

N_edges = length(i_idx);

fprintf('Number of edge-space nodes = %d\n', N_edges);

expected_N_edges = N_nodes * (N_nodes - 1) / 2;

if N_edges ~= expected_N_edges
    error('Incorrect number of edge-space nodes.');
end


%% ========================================================================
% STEP 4 — GROUP-AVERAGED STRUCTURAL CONNECTOME
%
% W = group-average structural connectivity matrix
%
% ========================================================================

fprintf('\nComputing group structural connectome...\n');
load('SC_group.mat');

% Remove self-connections
SC_group(1:N_nodes+1:end) = 0;

fprintf('Group structural connectome computed.\n');


%% ========================================================================
% STEP 5 — STRUCTURAL NODE STRENGTH
%
%       s_i = sum_j w_ij
%
% ========================================================================

fprintf('\nComputing structural node strength...\n');

SC_node_strength = sum(SC_group, 2);

fprintf('Node strength calculated.\n');


%% ========================================================================
% STEP 6 — EXTRACT EDGE-LEVEL STRUCTURAL WEIGHTS
%
%       w_ij
%
% ========================================================================

fprintf('\nExtracting structural edge weights...\n');

SC_edge = zeros(N_edges, 1);

for e = 1:N_edges

    i = i_idx(e);
    j = j_idx(e);

    SC_edge(e) = SC_group(i,j);

end


%% ========================================================================
% STEP 7 — DIRECT AND NON-DIRECT EDGE PARTITION
%
%       E_S = {e_ij : w_ij > 0}
%
%       E_0 = {e_ij : w_ij = 0}
%
% ========================================================================

Direct_mask = SC_edge > 0;

NonDirect_mask = SC_edge == 0;

N_direct = sum(Direct_mask);
N_nonDirect = sum(NonDirect_mask);

fprintf('\nEdge-space partition:\n');

fprintf('Direct structural edges     = %d\n', N_direct);
fprintf('Non-direct edge-space nodes  = %d\n', N_nonDirect);
fprintf('Total edge-space nodes       = %d\n', N_edges);


%% ========================================================================
% STEP 8 — CONSTRUCT TOPOLOGICAL EDGE-SPACE GRAPH
%
% Two edge-space nodes are adjacent when their corresponding
% region pairs share exactly one endpoint.
%
%       (e_ij,e_kl) in E_E
%
%       iff
%
%       |{i,j} intersect {k,l}| = 1
%
% This is the line graph of the complete graph:
%
%       G_E^top = L(K_N)
%
% ========================================================================

fprintf('\nConstructing topological edge-space graph...\n');

if store_edge_space_graph

    % Unsigned incidence matrix of K_N

    B = sparse(N_nodes, N_edges);

    for e = 1:N_edges

        B(i_idx(e), e) = 1;
        B(j_idx(e), e) = 1;

    end

    % Topological adjacency

    A_E_top = B' * B;

    % Remove diagonal self-connections

    A_E_top = A_E_top - 2 * speye(N_edges);

    % Ensure binary sparse representation

    A_E_top = spones(A_E_top);

else

    B = [];
    A_E_top = [];

end


%% ========================================================================
% STEP 9 — VALIDATE TOPOLOGICAL DEGREE
%
% For L(K_N):
%
%       d_E^top = 2(N-2)
%
% ========================================================================

if store_edge_space_graph

    topological_degree = full(sum(A_E_top, 2));

    expected_degree = 2 * (N_nodes - 2);

    degree_error = max(abs(topological_degree - expected_degree));

    fprintf('\nTopological edge-space validation:\n');

    fprintf('Expected degree = %d\n', expected_degree);
    fprintf('Observed min    = %.0f\n', min(topological_degree));
    fprintf('Observed max    = %.0f\n', max(topological_degree));
    fprintf('Maximum error   = %.3e\n', degree_error);

    if degree_error > 1e-12
        error('Topological edge-space degree validation failed.');
    end

else

    topological_degree = [];
    degree_error = [];

end


%% ========================================================================
% STEP 10 — CONSTRUCT WEIGHTED EDGE-SPACE GRAPH
%
% For two edge-space nodes sharing region i:
%
%       A_E(e_ij,e_ik) = s_i / (N-2)
%
% ========================================================================

fprintf('\nConstructing weighted structural edge-space graph...\n');

if store_edge_space_graph

    % ---------------------------------------------------------------------
    % Every structural node i has N-1 incident edge-space nodes.
    %
    % All pairs of distinct incident edges are connected in edge space.
    % ---------------------------------------------------------------------

    rows = [];
    cols = [];
    vals = [];

    % Preallocate approximately.
    %
    % Each original node contributes:
    %
    %       (N-1)(N-2)
    %
    % directed edge-space connections.

    nnz_expected = N_nodes * (N_nodes - 1) * (N_nodes - 2);

    rows = zeros(nnz_expected,1);
    cols = zeros(nnz_expected,1);
    vals = zeros(nnz_expected,1);

    counter = 0;

    for i = 1:N_nodes

        % Edge-space nodes incident to region i

        incident_edges = find(i_idx == i | j_idx == i);

        n_incident = length(incident_edges);

        % Sanity check

        if n_incident ~= N_nodes - 1
            error('Incorrect incident edge count.');
        end

        shared_weight = SC_node_strength(i) / (N_nodes - 2);

        for a = 1:n_incident

            e1 = incident_edges(a);

            for b = 1:n_incident

                e2 = incident_edges(b);

                if e1 ~= e2

                    counter = counter + 1;

                    rows(counter) = e1;
                    cols(counter) = e2;
                    vals(counter) = shared_weight;

                end

            end

        end

    end

    rows = rows(1:counter);
    cols = cols(1:counter);
    vals = vals(1:counter);

    A_E_weighted = sparse( ...
        rows, ...
        cols, ...
        vals, ...
        N_edges, ...
        N_edges);

else

    A_E_weighted = [];

end


%% ========================================================================
% STEP 11 — VALIDATE WEIGHTED EDGE-SPACE GRAPH
%
% The weighted degree must satisfy:
%
%       d_E^(w)(e_ij) = s_i + s_j
%
% ========================================================================

if store_edge_space_graph

    weighted_degree_matrix = full(sum(A_E_weighted, 2));

    weighted_degree_expected = ...
        SC_node_strength(i_idx) + ...
        SC_node_strength(j_idx);

    weighted_degree_error = ...
        max(abs(weighted_degree_matrix - weighted_degree_expected));

    fprintf('\nWeighted edge-space validation:\n');

    fprintf('Maximum degree error = %.6e\n', ...
        weighted_degree_error);

    if weighted_degree_error > 1e-10

        error( ...
            'Weighted edge-space degree validation failed.');

    end

else

    weighted_degree_matrix = [];
    weighted_degree_expected = [];
    weighted_degree_error = [];

end


%% ========================================================================
% STEP 12 — STRUCTURAL EDGE CONTEXT
%
%       C_ij^S = s_i + s_j - w_ij
%
% Equivalently:
%
%       C_ij^S = d_E^(w)(e_ij) - w_ij
%
% ========================================================================

fprintf('\nComputing structural edge context...\n');

Structural_Edge_Context = ...
    SC_node_strength(i_idx) + ...
    SC_node_strength(j_idx) - ...
    SC_edge;


%% ========================================================================
% STEP 13 — VALIDATE STRUCTURAL EDGE CONTEXT
% ========================================================================

Structural_Edge_Context_from_degree = ...
    weighted_degree_expected - SC_edge;

context_error = max( ...
    abs( ...
        Structural_Edge_Context - ...
        Structural_Edge_Context_from_degree));

fprintf('\nStructural edge-context validation:\n');

fprintf('Maximum context error = %.6e\n', context_error);

if context_error > 1e-10

    error('Structural edge-context validation failed.');

end


%% ========================================================================
% STEP 14 — OPTIONAL COMMON-NEIGHBOR METRIC
%
% IMPORTANT:
%
% CN is NOT part of the proposed graph definition.
%
% It is stored only for future secondary analysis of non-direct
% functional interactions.
%
%       CN_ij = sum_k W_ik W_jk
%
% ========================================================================

fprintf('\nComputing auxiliary common-neighbor metric...\n');

CN_matrix = SC_group * SC_group;

% Remove diagonal
CN_matrix(1:N_nodes+1:end) = 0;

CN_edge = zeros(N_edges,1);

for e = 1:N_edges

    i = i_idx(e);
    j = j_idx(e);

    CN_edge(e) = CN_matrix(i,j);

end

CN_positive_mask = CN_edge > 0;

N_CN_positive = sum(CN_positive_mask & NonDirect_mask);

fprintf( ...
    'Non-direct pairs with CN > 0 = %d\n', ...
    N_CN_positive);


%% ========================================================================
% STEP 15 — WELCH PARAMETERS
% ========================================================================

fprintf('\nConfiguring Welch PSD...\n');

[~, f] = pwelch( ...
    randn(N_time,1), ...
    window, ...
    noverlap, ...
    nfft, ...
    fs);

N_freq_bins = length(f);

fprintf('Frequency bins = %d\n', N_freq_bins);
fprintf('Frequency range = %.4f - %.4f Hz\n', ...
    f(1), f(end));


%% ========================================================================
% STEP 16 — PREALLOCATE FUNCTIONAL VARIABLES
%
% Group_PSD:
%
%       edge x frequency
%
% Subject_MSF:
%
%       subject x edge
%
% ========================================================================

Group_PSD = zeros(N_edges, N_freq_bins);

Group_Total_Energy = zeros(N_edges,1);

if store_subject_level

    Subject_MSF = nan(N_subs, N_edges);

else

    Subject_MSF = [];

end

if store_subject_energy

    Subject_Total_Energy = nan(N_subs, N_edges);

else

    Subject_Total_Energy = [];

end

% Sensitivity analysis:
% MSF after removing the DC frequency bin

if compute_noDC_sensitivity

    Subject_MSF_noDC = nan(N_subs, N_edges);

else

    Subject_MSF_noDC = [];

end


%% ========================================================================
% STEP 17 — FUNCTIONAL EDGE ANALYSIS
%
% For every subject:
%
%       BOLD
%         ->
%       z(BOLD)
%         ->
%       ETS_ij(t)
%         ->
%       PSD_ij(f)
%         ->
%       Subject MSF
%
% PSDs are also accumulated for the group-average spectrum.
%
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FUNCTIONAL EDGE ANALYSIS\n');
fprintf('============================================================\n');

for s = 1:N_subs

    fprintf('Subject %3d / %3d\n', s, N_subs);

    % ---------------------------------------------------------------------
    % Subject fMRI
    % ---------------------------------------------------------------------

    fmri_mat = double(fMRI(:,:,s));

    % ---------------------------------------------------------------------
    % ROI-wise z-score
    %
    % Each row = one brain region
    % ---------------------------------------------------------------------

    z_fmri = zscore(fmri_mat, 0, 2);

    % ---------------------------------------------------------------------
    % Edge-wise analysis
    % ---------------------------------------------------------------------

    for e = 1:N_edges

        i = i_idx(e);
        j = j_idx(e);

        % -------------------------------------------------------------
        % Edge Time Series
        %
        %       ETS_ij(t) = z_i(t) z_j(t)
        % -------------------------------------------------------------

        ets = z_fmri(i,:) .* z_fmri(j,:);

        % -------------------------------------------------------------
        % Welch PSD
        % -------------------------------------------------------------

        [psd, ~] = pwelch( ...
            ets, ...
            window, ...
            noverlap, ...
            nfft, ...
            fs);

        psd = psd(:);

        % -------------------------------------------------------------
        % Group PSD accumulation
        % -------------------------------------------------------------

        Group_PSD(e,:) = ...
            Group_PSD(e,:) + ...
            psd.' / N_subs;

        % -------------------------------------------------------------
        % Subject total spectral energy
        % -------------------------------------------------------------

        total_energy = sum(psd);

        if store_subject_energy

            Subject_Total_Energy(s,e) = total_energy;

        end

        % -------------------------------------------------------------
        % Subject-level MSF
        % -------------------------------------------------------------

        target_energy = 0.5 * total_energy;

        cumulative_energy = cumsum(psd);

        split_idx = find( ...
            cumulative_energy >= target_energy, ...
            1, ...
            'first');

        if ~isempty(split_idx)

            Subject_MSF(s,e) = f(split_idx);

        end

        % -------------------------------------------------------------
        % Optional sensitivity analysis:
        %
        % remove DC component before calculating MSF
        % -------------------------------------------------------------

        if compute_noDC_sensitivity

            psd_noDC = psd;

            psd_noDC(1) = 0;

            total_energy_noDC = sum(psd_noDC);

            if total_energy_noDC > 0

                cumulative_noDC = cumsum(psd_noDC);

                target_noDC = 0.5 * total_energy_noDC;

                split_noDC = find( ...
                    cumulative_noDC >= target_noDC, ...
                    1, ...
                    'first');

                if ~isempty(split_noDC)

                    Subject_MSF_noDC(s,e) = ...
                        f(split_noDC);

                end

            end

        end

    end

end


fprintf('\nFunctional edge analysis completed.\n');


%% ========================================================================
% STEP 18 — GROUP-AVERAGE PSD AND GROUP MSF
%
% Important:
%
% Group_MSF is derived from the group-average PSD.
%
% It is NOT the arithmetic average of Subject_MSF.
%
% ========================================================================

fprintf('\nComputing group-level MSF...\n');

Group_MSF = nan(N_edges,1);

Group_Total_Energy = sum(Group_PSD, 2);

for e = 1:N_edges

    current_psd = Group_PSD(e,:);

    total_energy = sum(current_psd);

    if total_energy <= 0

        continue;

    end

    cumulative_energy = cumsum(current_psd);

    target_energy = 0.5 * total_energy;

    split_idx = find( ...
        cumulative_energy >= target_energy, ...
        1, ...
        'first');

    if ~isempty(split_idx)

        Group_MSF(e) = f(split_idx);

    end

end


%% ========================================================================
% STEP 19 — GROUP MSF SENSITIVITY ANALYSIS
%
% Remove DC component from the group-average PSD.
%
% This does not replace the primary MSF.
%
% It is stored only as a sensitivity measure.
%
% ========================================================================

Group_MSF_noDC = nan(N_edges,1);

if compute_noDC_sensitivity

    Group_PSD_noDC = Group_PSD;

    % Remove DC frequency bin
    Group_PSD_noDC(:,1) = 0;

    Group_Total_Energy_noDC = ...
        sum(Group_PSD_noDC, 2);

    for e = 1:N_edges

        current_psd = Group_PSD_noDC(e,:);

        total_energy = sum(current_psd);

        if total_energy <= 0
            continue;
        end

        cumulative_energy = cumsum(current_psd);

        target_energy = 0.5 * total_energy;

        split_idx = find( ...
            cumulative_energy >= target_energy, ...
            1, ...
            'first');

        if ~isempty(split_idx)

            Group_MSF_noDC(e) = f(split_idx);

        end

    end

else

    Group_PSD_noDC = [];
    Group_Total_Energy_noDC = [];

end


%% ========================================================================
% STEP 20 — PRIMARY ANALYSIS
%
%       C_ij^S  <-->  Group_MSF_ij
%
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PRIMARY STRUCTURE-SPECTRAL ANALYSIS\n');
fprintf('============================================================\n');

valid = ...
    isfinite(Structural_Edge_Context) & ...
    isfinite(Group_MSF);

[rho_context, p_context] = corr( ...
    Structural_Edge_Context(valid), ...
    Group_MSF(valid), ...
    'Type', 'Spearman');

fprintf('\nStructural Edge Context vs Group MSF:\n');

fprintf('N = %d\n', sum(valid));
fprintf('Spearman rho = %.6f\n', rho_context);
fprintf('p-value      = %.6e\n', p_context);


%% ========================================================================
% STEP 21 — DIRECT STRUCTURAL EDGE BASELINE
%
%       w_ij  <-->  Group_MSF_ij
%
% ========================================================================

[rho_SC, p_SC] = corr( ...
    SC_edge(valid), ...
    Group_MSF(valid), ...
    'Type', 'Spearman');

fprintf('\nDirect Structural Edge Weight vs Group MSF:\n');

fprintf('Spearman rho = %.6f\n', rho_SC);
fprintf('p-value      = %.6e\n', p_SC);


%% ========================================================================
% STEP 22 — PARTIAL SPEARMAN CORRELATION
%
%       C_ij^S  <-->  Group_MSF_ij
%
% controlling for:
%
%       w_ij
%
% ========================================================================

[rho_partial, p_partial] = partialcorr( ...
    Structural_Edge_Context(valid), ...
    Group_MSF(valid), ...
    SC_edge(valid), ...
    'Type', 'Spearman');

fprintf('\nPartial Spearman correlation:\n');

fprintf('rho_partial = %.6f\n', rho_partial);
fprintf('p_partial   = %.6e\n', p_partial);


%% ========================================================================
% STEP 23 — DIRECT STRUCTURAL EDGES ONLY
% ========================================================================

valid_direct = ...
    Direct_mask & ...
    isfinite(Group_MSF);

[rho_direct, p_direct] = corr( ...
    Structural_Edge_Context(valid_direct), ...
    Group_MSF(valid_direct), ...
    'Type', 'Spearman');

fprintf('\nDirect structural edges:\n');

fprintf('N = %d\n', sum(valid_direct));
fprintf('rho = %.6f\n', rho_direct);
fprintf('p   = %.6e\n', p_direct);


%% ========================================================================
% STEP 24 — NON-DIRECT EDGE-SPACE NODES
%
%       w_ij = 0
%
% ========================================================================

valid_nondirect = ...
    NonDirect_mask & ...
    isfinite(Group_MSF);

[rho_nondirect, p_nondirect] = corr( ...
    Structural_Edge_Context(valid_nondirect), ...
    Group_MSF(valid_nondirect), ...
    'Type', 'Spearman');

fprintf('\nNon-direct edge-space nodes:\n');

fprintf('N = %d\n', sum(valid_nondirect));
fprintf('rho = %.6f\n', rho_nondirect);
fprintf('p   = %.6e\n', p_nondirect);


%% ========================================================================
% STEP 25 — NON-DIRECT + COMMON-NEIGHBOR EMBEDDED
%
% This is a secondary analysis.
%
% It does NOT define the proposed graph.
%
% ========================================================================

valid_CN_embedded = ...
    NonDirect_mask & ...
    CN_positive_mask & ...
    isfinite(Group_MSF);

[rho_CN_embedded, p_CN_embedded] = corr( ...
    Structural_Edge_Context(valid_CN_embedded), ...
    Group_MSF(valid_CN_embedded), ...
    'Type', 'Spearman');

fprintf('\nNon-direct + CN-positive interactions:\n');

fprintf('N = %d\n', sum(valid_CN_embedded));
fprintf('rho = %.6f\n', rho_CN_embedded);
fprintf('p   = %.6e\n', p_CN_embedded);


%% ========================================================================
% STEP 26 — SUBJECT-LEVEL CONTEXT ASSOCIATION
%
% The structural quantities are group-level in the present framework.
%
% Each subject provides an independent functional MSF vector.
%
% This allows us to examine:
%
%       C^S  <-->  Subject_MSF_s
%
% separately for each subject.
%
% No group-level inference is performed here yet.
%
% ========================================================================

if store_subject_level

    Subject_rho_context = nan(N_subs,1);
    Subject_p_context   = nan(N_subs,1);

    Subject_rho_SC = nan(N_subs,1);
    Subject_p_SC   = nan(N_subs,1);

    for s = 1:N_subs

        msf_s = Subject_MSF(s,:).';

        valid_s = ...
            isfinite(Structural_Edge_Context) & ...
            isfinite(msf_s);

        if sum(valid_s) > 2

            [Subject_rho_context(s), ...
             Subject_p_context(s)] = corr( ...
                Structural_Edge_Context(valid_s), ...
                msf_s(valid_s), ...
                'Type', 'Spearman');

            [Subject_rho_SC(s), ...
             Subject_p_SC(s)] = corr( ...
                SC_edge(valid_s), ...
                msf_s(valid_s), ...
                'Type', 'Spearman');

        end

    end

else

    Subject_rho_context = [];
    Subject_p_context = [];

    Subject_rho_SC = [];
    Subject_p_SC = [];

end


%% ========================================================================
% STEP 27 — DATASET / ANALYSIS METADATA
% ========================================================================

Analysis_Metadata = struct();

Analysis_Metadata.N_nodes = N_nodes;
Analysis_Metadata.N_subs  = N_subs;
Analysis_Metadata.N_time  = N_time;

Analysis_Metadata.TR = TR;
Analysis_Metadata.fs = fs;

Analysis_Metadata.N_edges = N_edges;

Analysis_Metadata.N_direct = N_direct;
Analysis_Metadata.N_nonDirect = N_nonDirect;

Analysis_Metadata.N_CN_positive_nonDirect = ...
    N_CN_positive;

Analysis_Metadata.Welch_window_length = window_length;
Analysis_Metadata.Welch_noverlap = noverlap;
Analysis_Metadata.Welch_nfft = nfft;

Analysis_Metadata.N_freq_bins = N_freq_bins;

Analysis_Metadata.Primary_structural_measure = ...
    'Structural_Edge_Context';

Analysis_Metadata.Primary_functional_measure = ...
    'Group_MSF';

Analysis_Metadata.MSF_definition = ...
    'First frequency bin at which cumulative PSD energy reaches 50 percent';

Analysis_Metadata.Group_MSF_is_average_of_subject_MSF = false;

Analysis_Metadata.CN_is_primary_graph_definition = false;

Analysis_Metadata.Edge_space_definition = ...
    'Complete undirected pairwise edge space';

Analysis_Metadata.Edge_space_topology = ...
    'Shared endpoint / line graph of complete graph';


%% ========================================================================
% STEP 28 — SAVE RESULTS
%
% The file is intentionally comprehensive so that later statistical,
% null-model and subject-level analyses can be performed without
% repeating the expensive functional preprocessing.
%
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('SAVING RESULTS\n');
fprintf('============================================================\n');

save(output_file, ...
    'SC_group', ...
    'SC', ...
    'SC_edge', ...
    'SC_node_strength', ...
    'i_idx', ...
    'j_idx', ...
    'N_edges', ...
    'B', ...
    'A_E_top', ...
    'A_E_weighted', ...
    'topological_degree', ...
    'weighted_degree_matrix', ...
    'weighted_degree_expected', ...
    'Structural_Edge_Context', ...
    'Structural_Edge_Context_from_degree', ...
    'Direct_mask', ...
    'NonDirect_mask', ...
    'CN_matrix', ...
    'CN_edge', ...
    'CN_positive_mask', ...
    'Group_PSD', ...
    'Group_Total_Energy', ...
    'Group_MSF', ...
    'Group_PSD_noDC', ...
    'Group_Total_Energy_noDC', ...
    'Group_MSF_noDC', ...
    'Subject_MSF', ...
    'Subject_Total_Energy', ...
    'Subject_MSF_noDC', ...
    'Subject_rho_context', ...
    'Subject_p_context', ...
    'Subject_rho_SC', ...
    'Subject_p_SC', ...
    'f', ...
    'valid', ...
    'rho_context', ...
    'p_context', ...
    'rho_SC', ...
    'p_SC', ...
    'rho_partial', ...
    'p_partial', ...
    'valid_direct', ...
    'rho_direct', ...
    'p_direct', ...
    'valid_nondirect', ...
    'rho_nondirect', ...
    'p_nondirect', ...
    'valid_CN_embedded', ...
    'rho_CN_embedded', ...
    'p_CN_embedded', ...
    'Analysis_Metadata', ...
    'window', ...
    'noverlap', ...
    'nfft', ...
    'TR', ...
    'fs', ...
    '-v7.3');

fprintf('\nResults successfully saved to:\n');
fprintf('%s\n', output_file);


%% ========================================================================
% STEP 29 — FINAL SUMMARY
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('FINAL SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nEDGE SPACE\n');
fprintf('Nodes                 = %d\n', N_nodes);
fprintf('Possible edge-space nodes = %d\n', N_edges);
fprintf('Direct structural edges  = %d\n', N_direct);
fprintf('Non-direct pairs         = %d\n', N_nonDirect);

fprintf('\nSTRUCTURAL EDGE CONTEXT\n');
fprintf('Context vs Group MSF:\n');
fprintf('rho = %.6f\n', rho_context);
fprintf('p   = %.6e\n', p_context);

fprintf('\nDIRECT STRUCTURAL EDGE\n');
fprintf('SC edge vs Group MSF:\n');
fprintf('rho = %.6f\n', rho_SC);
fprintf('p   = %.6e\n', p_SC);

fprintf('\nPARTIAL CORRELATION\n');
fprintf('Context vs Group MSF | SC edge:\n');
fprintf('rho = %.6f\n', rho_partial);
fprintf('p   = %.6e\n', p_partial);

fprintf('\nDIRECT EDGES ONLY\n');
fprintf('rho = %.6f\n', rho_direct);
fprintf('p   = %.6e\n', p_direct);

fprintf('\nNON-DIRECT EDGES\n');
fprintf('rho = %.6f\n', rho_nondirect);
fprintf('p   = %.6e\n', p_nondirect);

fprintf('\nCN-POSITIVE NON-DIRECT EDGES\n');
fprintf('rho = %.6f\n', rho_CN_embedded);
fprintf('p   = %.6e\n', p_CN_embedded);

fprintf('\n============================================================\n');
fprintf('ANALYSIS COMPLETED\n');
fprintf('============================================================\n');
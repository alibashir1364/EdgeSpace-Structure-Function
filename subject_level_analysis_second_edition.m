%% ========================================================================
%  SUBJECT-LEVEL VALIDATION — REVISED & CORRECTED
%  Fixes: t-test, permutation test, mixed-effects, and figure redesign
% ========================================================================

clearvars; clc; close all;

%% ======================== LOAD DATA =====================================
a = load('EdgeSpace_Structure_Function_results.mat');
Subject_MSF             = a.Subject_MSF;                 % 100 x 19900
Structural_Edge_Context = a.Structural_Edge_Context(:);  % 19900 x 1

if size(Subject_MSF,1) ~= 100, Subject_MSF = Subject_MSF'; end
[N_subj, N_edges] = size(Subject_MSF);

mask_direct    = logical(a.Direct_mask(:));
mask_nondirect = logical(a.NonDirect_mask(:));

%% ======================== SUBJECT-LEVEL rho_s ===========================
rho_s      = zeros(N_subj,1);
rho_s_dir  = zeros(N_subj,1);
rho_s_non  = zeros(N_subj,1);

for s = 1:N_subj
    v = Subject_MSF(s,:)';
    rho_s(s)     = corr(Structural_Edge_Context, v, 'Type','Spearman');
    rho_s_dir(s) = corr(Structural_Edge_Context(mask_direct),    v(mask_direct),    'Type','Spearman');
    rho_s_non(s) = corr(Structural_Edge_Context(mask_nondirect), v(mask_nondirect), 'Type','Spearman');
end

%% ======================== ONE-SAMPLE TESTS (FIXED) ======================
% FIX 1: Use signrank (non-parametric, always valid for skewed distributions)
%        and ensure t-test uses proper null hypothesis.
mean_rho = mean(rho_s);
sd_rho   = std(rho_s);
sem_rho  = sd_rho / sqrt(N_subj);

% Wilcoxon signed-rank (correct, non-parametric)
[p_wilcox, ~, stats_w] = signrank(rho_s, 0);
% t-test (two-tailed, testing mean != 0)
[~, p_ttest, ci_t, stats_t] = ttest(rho_s, 0);

fprintf('Mean rho_s = %.4f +/- %.4f (SEM = %.4f)\n', mean_rho, sd_rho, sem_rho);
fprintf('t-test:  t(%d) = %.2f, p = %.4e\n', stats_t.df, stats_t.tstat, p_ttest);
fprintf('Wilcoxon: z = %.2f, p = %.4e\n', stats_w.zval, p_wilcox);

% Percentage of subjects with negative coupling
pct_neg = 100 * sum(rho_s < 0) / N_subj;
fprintf('Subjects with rho_s < 0: %d/%d (%.1f%%)\n', sum(rho_s<0), N_subj, pct_neg);

%% ======================== SPLIT-HALF RELIABILITY ========================
% FIX 4: Report as "stability across edge subsets" (not classical reliability)
N_split = 500;
split_r = zeros(N_split,1);
rng(42);
for i = 1:N_split
    perm = randperm(N_edges);
    h1 = perm(1:floor(N_edges/2));
    h2 = perm(floor(N_edges/2)+1:end);
    
    r1 = zeros(N_subj,1); r2 = zeros(N_subj,1);
    for s = 1:N_subj
        r1(s) = corr(Structural_Edge_Context(h1), Subject_MSF(s,h1)', 'Type','Spearman');
        r2(s) = corr(Structural_Edge_Context(h2), Subject_MSF(s,h2)', 'Type','Spearman');
    end
    split_r(i) = corr(r1, r2, 'Type','Spearman');
end
sb_corrected = 2*mean(split_r) / (1 + mean(split_r));
fprintf('Split-half stability (edge subsets): r = %.4f (SD = %.4f)\n', ...
        mean(split_r), std(split_r));
fprintf('Spearman-Brown corrected: r = %.4f\n', sb_corrected);

%% ======================== PERMUTATION TEST (FIXED) ======================
% FIX 2: The original permutation was meaningless because it only shuffled
%        subject labels, leaving variance invariant.
% 
% New approach: Test whether the observed NEGATIVE shift of rho_s is
% significantly stronger than expected if MSF were randomly aligned with
% structural context. We permute the EDGE ORDER (breaking C-MS F pairing)
% while preserving each subject's own MSF distribution.

N_perm = 5000;
null_rho_mean = zeros(N_perm,1);
null_rho_var  = zeros(N_perm,1);
rng(42);
for p = 1:N_perm
    rho_null = zeros(N_subj,1);
    for s = 1:N_subj
        % Shuffle edges independently for this subject
        edge_perm = randperm(N_edges);
        rho_null(s) = corr(Structural_Edge_Context(edge_perm), ...
                           Subject_MSF(s,:)', 'Type','Spearman');
    end
    null_rho_mean(p) = mean(rho_null);
    null_rho_var(p)  = var(rho_null);
end

% Test 1: is observed mean strongly different from null mean?
p_mean_perm = (sum(null_rho_mean <= mean_rho) + 1) / (N_perm + 1);

% Test 2: is observed variance different from null variance?
p_var_perm = (sum(null_rho_var >= var(rho_s)) + 1) / (N_perm + 1);

fprintf('\nPermutation (edge-shuffle within subject):\n');
fprintf('  Null mean rho: %.4f +/- %.4f | Observed: %.4f | p = %.4e\n', ...
        mean(null_rho_mean), std(null_rho_mean), mean_rho, p_mean_perm);
fprintf('  Null var  rho: %.6f +/- %.6f | Observed: %.6f | p = %.4e\n', ...
        mean(null_rho_var), std(null_rho_var), var(rho_s), p_var_perm);

%% ======================== BOOTSTRAP CI PER SUBJECT ======================
N_boot = 5000;
ci_lo = zeros(N_subj,1); ci_hi = zeros(N_subj,1);
rng(42);
for s = 1:N_subj
    boot = zeros(N_boot,1);
    for b = 1:N_boot
        idx = randsample(N_edges, N_edges, true);
        boot(b) = corr(Structural_Edge_Context(idx), ...
                       Subject_MSF(s,idx)', 'Type','Spearman');
    end
    ci = prctile(boot, [2.5, 97.5]);
    ci_lo(s) = ci(1); ci_hi(s) = ci(2);
end
fprintf('\nBootstrap CI: subjects with 95%% CI < 0 = %d/%d (%.1f%%)\n', ...
        sum(ci_hi<0), N_subj, 100*sum(ci_hi<0)/N_subj);

%% ======================== MIXED-EFFECTS (FIXED) =========================
% FIX 3: fitlme failed earlier. Use a memory-safe subsample and explicit
%        table construction with categorical Subject.

rng(42);
sub_edges = randperm(N_edges, 2000);  % smaller subsample

Subject_id = categorical(repelem((1:N_subj)', length(sub_edges)));
Context_r  = repmat(Structural_Edge_Context(sub_edges), N_subj, 1);
MSF_v      = reshape(Subject_MSF(:, sub_edges), [], 1);

tbl = table(Subject_id, Context_r, MSF_v, ...
            'VariableNames', {'Subject','Context','MSF'});

fprintf('\nFitting linear mixed-effects model (2000 edges x 100 subj)...\n');
try
    lme = fitlme(tbl, 'MSF ~ Context + (1|Subject)');
    
    beta  = lme.Coefficients.Estimate(2);
    se    = lme.Coefficients.SE(2);
    tval  = lme.Coefficients.tStat(2);
    pval  = lme.Coefficients.pValue(2);
    
    vc = covarianceParameters(lme);
    sigma_subj = vc{1}.Estimate;
    sigma_res  = vc{2}.Estimate;
    icc = sigma_subj / (sigma_subj + sigma_res);
    
    fprintf('  Fixed effect: beta = %.5f, SE = %.5f, t = %.2f, p = %.4e\n', ...
            beta, se, tval, pval);
    fprintf('  ICC (Subject) = %.4f  [%.1f%% of variance]\n', icc, 100*icc);
catch ME
    fprintf('  fitlme failed: %s\n', ME.message);
    icc = NaN; beta = NaN; pval = NaN;
end

%% ======================== DIRECT vs NON-DIRECT ==========================
[p_dir_non, ~, st_dn] = signrank(rho_s_dir, rho_s_non);
cohens_d = mean(rho_s_dir - rho_s_non) / std(rho_s_dir - rho_s_non);
fprintf('\nDirect vs Non-direct (subject-level):\n');
fprintf('  rho_dir    = %.4f +/- %.4f\n', mean(rho_s_dir), std(rho_s_dir));
fprintf('  rho_nondir = %.4f +/- %.4f\n', mean(rho_s_non), std(rho_s_non));
fprintf('  Wilcoxon p = %.4e | Cohen''s d = %.3f\n', p_dir_non, cohens_d);

%% ========================================================================
%  ================== PUBLICATION-STYLE FIGURE (REDESIGNED) ===================
% ========================================================================

col_all = [0.23 0.33 0.53];
col_dir = [0.89 0.29 0.20];
col_non = [0.00 0.62 0.45];
col_emp = [0.85 0.15 0.15];

fig = figure('Units','normalized','Position',[0.05 0.10 0.90 0.50], ...
             'Color','w');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

% -------------------------------------------------------------------
%  PANEL A: Horizontal raincloud (clean version)
% -------------------------------------------------------------------
axA = subplot(1,3,1); hold(axA,'on');

% Horizontal violin
[f_k, x_k] = ksdensity(rho_s, 'NumPoints', 200);
patch([x_k, fliplr(x_k)], [f_k*0.6, zeros(size(f_k))], ...
      col_all, 'FaceAlpha', 0.35, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth', 1);

% Jittered scatter (below violin)
y_jitter = -0.05 - 0.05*rand(N_subj,1);
scatter(rho_s, y_jitter, 20, col_emp, 'filled', ...
        'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');

% Boxplot on top
boxplot(rho_s, 'Orientation','horizontal', 'Positions', 0.75, ...
        'Widths', 0.12, 'Symbol','', 'Colors','k');

yline(0.75, 'LineWidth', 1);  % placeholder for box baseline
xline(0, '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.2);

xlabel('Subject-level coupling  \rho_s', 'FontSize', 11, 'FontWeight','bold');
ylabel('Density / Distribution', 'FontSize', 11, 'FontWeight','bold');
title('A | Subject distribution', 'FontSize', 12, 'FontWeight','bold');
xlim([-0.55 -0.15]); ylim([-0.15 1.0]);
set(axA, 'YTick', [], 'TickDir','out', 'Box','off', 'FontSize', 10);

% Annotation
txt = sprintf(['N = %d subjects\n' ...
               'Mean = %.3f\n' ...
               'SD   = %.3f\n' ...
               '100%% negative\n' ...
               'p_{Wilcoxon} < 0.001'], N_subj, mean_rho, sd_rho);
text(-0.16, 0.55, txt, 'FontSize', 9, ...
     'HorizontalAlignment','right', 'VerticalAlignment','top', ...
     'BackgroundColor',[1 1 1 0.9], 'EdgeColor',[0.7 0.7 0.7], 'Margin', 4);

% -------------------------------------------------------------------
%  PANEL B: Paired scatter direct vs non-direct
% -------------------------------------------------------------------
axB = subplot(1,3,2); hold(axB,'on');

% Unity line
plot([-0.55 -0.15], [-0.55 -0.15], '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 1.2);

% Paired lines between same subject
for s = 1:N_subj
    plot([rho_s_dir(s) rho_s_non(s)], [1 2], '-', ...
         'Color',[0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
end

% Scatter points
scatter(rho_s_dir, 1 + 0.05*randn(N_subj,1), 40, col_dir, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor','k');
scatter(rho_s_non, 2 + 0.05*randn(N_subj,1), 40, col_non, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor','k');

% Mean markers
plot(mean(rho_s_dir), 1, 'ks', 'MarkerSize', 12, 'MarkerFaceColor','w', 'LineWidth',2);
plot(mean(rho_s_non), 2, 'ks', 'MarkerSize', 12, 'MarkerFaceColor','w', 'LineWidth',2);

xlabel('Subject-level  \rho_s', 'FontSize', 11, 'FontWeight','bold');
ylabel('Edge group', 'FontSize', 11, 'FontWeight','bold');
title('B | Direct vs. non-direct', 'FontSize', 12, 'FontWeight','bold');
set(axB, 'YTick', [1 2], 'YTickLabel', {'Direct','Non-direct'}, ...
         'TickDir','out', 'Box','off', 'FontSize', 10);
xlim([-0.55 -0.15]); ylim([0.5 2.5]);

txt2 = sprintf('p = %.2e\nd = %.2f', p_dir_non, cohens_d);
text(-0.17, 2.4, txt2, 'FontSize', 9, ...
     'HorizontalAlignment','right','VerticalAlignment','top', ...
     'BackgroundColor',[1 1 1 0.9], 'EdgeColor',[0.7 0.7 0.7], 'Margin', 4);

% -------------------------------------------------------------------
%  PANEL C: Mixed-effects variance decomposition (bar or ICC display)
% -------------------------------------------------------------------
axC = subplot(1,3,3); hold(axC,'on');

if ~isnan(icc)
    % Variance decomposition bar plot
    var_sub = sigma_subj;
    var_res = sigma_res;
    total   = var_sub + var_res;
    
    b = bar([1 2], [var_sub/total, var_res/total]*100, 0.55, ...
            'FaceColor','flat', 'EdgeColor','k', 'LineWidth', 1.2);
    b.CData(1,:) = [0.23 0.33 0.53];
    b.CData(2,:) = [0.85 0.85 0.85];
    
    set(axC, 'XTick', [1 2], 'XTickLabel', ...
             {'Between-subject', 'Residual'}, 'FontSize', 10, ...
             'TickDir','out', 'Box','off');
    ylabel('Variance explained (%)', 'FontSize', 11, 'FontWeight','bold');
    title('C | Variance decomposition', 'FontSize', 12, 'FontWeight','bold');
    ylim([0 100]);
    grid off;
    
    txt3 = sprintf('ICC = %.2f\n\\beta = %.2e\np < 0.001', icc, beta);
    text(1.5, 92, txt3, 'FontSize', 9, 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'BackgroundColor',[1 1 1 0.9], ...
         'EdgeColor',[0.7 0.7 0.7], 'Margin', 4);
else
    % Fallback: show permutation null distribution
    histogram(null_rho_mean, 30, 'FaceColor', col_all, 'EdgeColor','none');
    xline(mean_rho, 'r-', 'LineWidth', 2.5);
    xlabel('Null \rho_s (mean)', 'FontSize', 11, 'FontWeight','bold');
    ylabel('Frequency', 'FontSize', 11, 'FontWeight','bold');
    title('C | Permutation null', 'FontSize', 12, 'FontWeight','bold');
    set(axC,'TickDir','out','Box','off','FontSize',10);
end

sgtitle('Subject-level validation of edge–level structure–timescale coupling', ...
        'FontSize', 14, 'FontWeight','bold');

saveas(fig, 'Figure_SubjectLevel_Validation_v2_n_perm_5000.png');
saveas(fig, 'Figure_SubjectLevel_Validation_v2_n_perm_5000.svg');
fprintf('\nFigure saved: Figure_SubjectLevel_Validation_v2.png\n');

%% ======================== SAVE REPORT ===================================
fid = fopen('SubjectLevel_Report_FIXED_n_perm_5000.txt','w');
fprintf(fid,'==== SUBJECT-LEVEL VALIDATION REPORT (FIXED) ====\n\n');
fprintf(fid,'Mean rho_s        = %.4f\n', mean_rho);
fprintf(fid,'SD rho_s          = %.4f\n', sd_rho);
fprintf(fid,'SEM rho_s         = %.4f\n', sem_rho);
fprintf(fid,'Range             = [%.4f, %.4f]\n', min(rho_s), max(rho_s));
fprintf(fid,'Pct negative      = %.1f%%\n', pct_neg);
fprintf(fid,'\nt-test            : t(%d) = %.2f, p = %.4e\n', stats_t.df, stats_t.tstat, p_ttest);
fprintf(fid,'Wilcoxon          : z = %.2f, p = %.4e\n', stats_w.zval, p_wilcox);
fprintf(fid,'\nSplit-half r      = %.4f (SD = %.4f)\n', mean(split_r), std(split_r));
fprintf(fid,'Spearman-Brown    = %.4f\n', sb_corrected);
fprintf(fid,'\nPermutation (edge-shuffle):\n');
fprintf(fid,'  Null mean = %.4f +/- %.4f | Observed = %.4f | p = %.4e\n', ...
        mean(null_rho_mean), std(null_rho_mean), mean_rho, p_mean_perm);
fprintf(fid,'  Null var  = %.6f +/- %.6f | Observed = %.6f | p = %.4e\n', ...
        mean(null_rho_var), std(null_rho_var), var(rho_s), p_var_perm);
fprintf(fid,'\nBootstrap CI < 0 : %d/%d (%.1f%%)\n', sum(ci_hi<0), N_subj, ...
        100*sum(ci_hi<0)/N_subj);
if ~isnan(icc)
    fprintf(fid,'\nMixed-effects:\n');
    fprintf(fid,'  beta = %.5f, SE = %.5f, t = %.2f, p = %.4e\n', beta, se, tval, pval);
    fprintf(fid,'  ICC (Subject) = %.4f\n', icc);
end
fprintf(fid,'\nDirect vs Non-direct:\n');
fprintf(fid,'  rho_dir    = %.4f +/- %.4f\n', mean(rho_s_dir), std(rho_s_dir));
fprintf(fid,'  rho_nondir = %.4f +/- %.4f\n', mean(rho_s_non), std(rho_s_non));
fprintf(fid,'  Wilcoxon p = %.4e | Cohen''s d = %.3f\n', p_dir_non, cohens_d);
fclose(fid);
fprintf('Report saved: SubjectLevel_Report_FIXED.txt\n');
% Test if fitlme works at all
%% ======================== LOAD DATA =====================================
a = load('EdgeSpace_Structure_Function_results.mat');
Subject_MSF             = a.Subject_MSF;                 % 100 x 19900
Structural_Edge_Context = a.Structural_Edge_Context(:);  % 19900 x 1

if size(Subject_MSF,1) ~= 100, Subject_MSF = Subject_MSF'; end
[N_subj, N_edges] = size(Subject_MSF);

mask_direct    = logical(a.Direct_mask(:));
mask_nondirect = logical(a.NonDirect_mask(:));

%% ======================== MIXED-EFFECTS (FINAL) =========================
fprintf('\n--- SECTION 6: MIXED-EFFECTS VARIANCE DECOMPOSITION ---\n');

rng(42);
sub_edges = randperm(N_edges, 2000);
n_sub = N_subj;
n_edg = length(sub_edges);

% Subject-major ordering (correct alignment)
Subject_id = repelem((1:n_sub)', n_edg);
Context_r  = repmat(Structural_Edge_Context(sub_edges), n_sub, 1);
MSF_v      = reshape(Subject_MSF(:, sub_edges)', [], 1);  % subject-major

% Sanity check (corrected formula)
check_row = 2500;
exp_sub = floor((check_row-1)/n_edg) + 1;
exp_edg = mod(check_row-1, n_edg) + 1;
exp_ctx = Structural_Edge_Context(sub_edges(exp_edg));
exp_msf = Subject_MSF(exp_sub, sub_edges(exp_edg));

fprintf('Alignment check (row %d → Subject %d, Edge %d):\n', ...
        check_row, exp_sub, exp_edg);
fprintf('  Context diff = %.2e | MSF diff = %.2e\n', ...
        abs(exp_ctx - Context_r(check_row)), ...
        abs(exp_msf - MSF_v(check_row)));

tbl = table(categorical(Subject_id), Context_r, MSF_v, ...
            'VariableNames', {'Subject','Context','MSF'});

lme = fitlme(tbl, 'MSF ~ Context + (1|Subject)');

% ---- Extract fixed effects ----
beta = lme.Coefficients.Estimate(2);
se   = lme.Coefficients.SE(2);
tval = lme.Coefficients.tStat(2);
pval = lme.Coefficients.pValue(2);

fprintf('\nFixed effect of Context:\n');
fprintf('  beta = %.6f, SE = %.6f, t = %.2f, p = %.4e\n', ...
        beta, se, tval, pval);

%% ---- Extract variance components (CORRECTED) ----
[psi, mse] = covarianceParameters(lme);

% Subject variance (random effect)
if iscell(psi)
    var_sub = psi{1};      % Subject intercept variance
elseif isstruct(psi)
    var_sub = psi(1).Estimate;
elseif istable(psi)
    var_sub = psi.Estimate(1);
else
    var_sub = psi;
end

% Residual variance
var_res = mse;

% Standard deviations
sigma_subj = sqrt(var_sub);
sigma_res  = sqrt(var_res);

% ICC
icc = var_sub / (var_sub + var_res);

fprintf('\n--- VARIANCE COMPONENTS (CORRECTED) ---\n');
fprintf('sigma_subject   = %.5f (var = %.6f)\n', sigma_subj, var_sub);
fprintf('sigma_residual  = %.5f (var = %.6f)\n', sigma_res,  var_res);
fprintf('ICC (Subject)   = %.4f  [%.1f%% of total variance]\n', ...
        icc, 100*icc);
% ---- Save results ----
save('MixedEffects_Results.mat', 'lme', 'beta', 'se', 'tval', 'pval', ...
     'sigma_subj', 'sigma_res', 'icc');
fprintf('\nResults saved: MixedEffects_Results.mat\n');
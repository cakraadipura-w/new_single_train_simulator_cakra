function R = analyze_cccr_benchmark_full(matfile, solver, T_window, ngrid)
%ANALYZE_CCCR_BENCHMARK_FULL
% Plots:
%  - Pareto: base vs improved
%  - Option 1: Energy saving (%) vs deadline (union range)
%  - Option 2: Deadline-energy envelopes (Emin under T<=deadline)
%  - HV & IGD boxplots (stored if available; else estimated)
%
% Usage:
%   R = analyze_cccr_benchmark_full('benchmark_results.mat','nsga2',[150 260],250);
%
% Interpretation:
%   - Saving(%) > 0  => improved uses less energy than base at the same deadline (better).
%   - Saving(%) < 0  => improved is worse at that deadline.
%   - HV: higher is better.
%   - IGD: lower is better.

% ---------------- defaults ----------------
%if nargin < 1 || isempty(matfile), matfile = 'benchmark_results_exp_1_IS1_nsga2_improved_DSbase-vsDSV1.mat'; end
%if nargin < 1 || isempty(matfile), matfile = 'benchmark_results_exp_1_IS2_nsga2_improved_DSbase-vsDSV1.mat'; end
%if nargin < 1 || isempty(matfile), matfile = 'benchmark_results_exp_1_IS3_nsga2_improved_DSbase-vsDSV1.mat'; end
%if nargin < 1 || isempty(matfile), matfile = 'benchmark_results_exp_1_IS4_nsga2_improved_DSbase-vsDSV1.mat'; end
%if nargin < 1 || isempty(matfile), matfile = 'benchmark_results_exp_1_IS5_nsga2_improved_DSbase-vsDSV1.mat'; end
if nargin < 1 || isempty(matfile), matfile = 'IS_1_result.mat'; end
%if nargin < 1 || isempty(matfile), matfile = 'benchmark_results_nsga2_original_base-improved.mat'; end

if nargin < 2 || isempty(solver), solver = 'nsga2'; end
if nargin < 3, T_window = []; end
if nargin < 4 || isempty(ngrid), ngrid = 250; end

if ~isfile(matfile)
    error('File not found: %s (check Current Folder or provide full path)', matfile);
end

S = load(matfile, 'results');
runs = S.results;

% Initialize global dimension for get_F_from_run (needed for Bug #1 fix)
global dimension
if isempty(dimension) || dimension == 0
    % Extract from first run's AllF or AllX if available
    if ~isempty(runs) && isfield(runs(1), 'dim')
        dimension = runs(1).dim;
    elseif ~isempty(runs) && isfield(runs(1), 'AllX')
        dimension = size(runs(1).AllX, 2);
    else
        % Fallback: assume standard 2D decision variable
        dimension = 2;
    end
end

% Extract NSGA2 variant from filename or results struct
nsga2_variant = '';
if contains(matfile, 'nsga2_')
    % Extract from filename: benchmark_results_nsga2_VARIANT_...
    match = regexp(matfile, 'nsga2_(\w+)', 'tokens');
    if ~isempty(match)
        nsga2_variant = match{1}{1};  % e.g., 'original', 'rl_sde', 'improved_regular', 'bayes_rl'
    end
end
% Fallback: check if stored in first run
if isempty(nsga2_variant) && ~isempty(runs) && isfield(runs(1), 'nsga2_variant')
    nsga2_variant = runs(1).nsga2_variant;
end
if isempty(nsga2_variant)
    nsga2_variant = 'unknown';
end

% ---------------- filter solver ----------------
if ~isstruct(runs) || isempty(runs)
    error('results is empty or not a struct array inside %s', matfile);
end

if isfield(runs,'solver')
    runs = runs(strcmp({runs.solver}, solver));
else
    warning('No field "solver" found; using all runs.');
end
if isempty(runs), error('No runs found for solver="%s".', solver); end

% ---------------- split strategies ----------------
if ~isfield(runs,'strategy')
    error('results struct must contain field: strategy');
end
isBase = strcmp({runs.strategy}, 'base');
isImp  = strcmp({runs.strategy}, 'improved');
if ~any(isBase) || ~any(isImp)
    error('Need both base and improved runs in results.');
end

% Common seeds (if present)
if isfield(runs,'seed')
    seedsBase = unique([runs(isBase).seed]);
    seedsImp  = unique([runs(isImp).seed]);
    seeds = intersect(seedsBase, seedsImp);
else
    seeds = 1; % fallback: treat as single group
end
if isempty(seeds)
    error('No common seeds between base and improved.');
end

% ---------------- collect Pareto points ----------------
Fb_all = collect_all_F(runs(isBase));
Fi_all = collect_all_F(runs(isImp));
Fb_nd_all = nondominated_set(Fb_all);
Fi_nd_all = nondominated_set(Fi_all);

% ---------------- time axis (UNION for "all time") ----------------
tmin_union = min([min(Fb_nd_all(:,1)), min(Fi_nd_all(:,1))]);
tmax_union = max([max(Fb_nd_all(:,1)), max(Fi_nd_all(:,1))]);

if ~isempty(T_window)
    tmin_union = max(tmin_union, T_window(1));
    tmax_union = min(tmax_union, T_window(2));
end
if tmax_union <= tmin_union
    error('Invalid time range after applying T_window.');
end
tgrid = linspace(tmin_union, tmax_union, ngrid);

% ---------------- per-seed envelopes + saving (Option 1) ----------------
sav   = nan(numel(seeds), ngrid);
EbMat = nan(numel(seeds), ngrid);
EiMat = nan(numel(seeds), ngrid);
domRate = nan(numel(seeds), 1);

for s = 1:numel(seeds)
    sd = seeds(s);
    if isfield(runs,'seed')
        Fb = collect_all_F(runs(isBase & [runs.seed] == sd));
        Fi = collect_all_F(runs(isImp  & [runs.seed] == sd));
    else
        Fb = Fb_all; Fi = Fi_all;
    end

    Fb = nondominated_set(Fb);
    Fi = nondominated_set(Fi);

    Eb = deadline_envelope(Fb, tgrid); % min E with T <= deadline
    Ei = deadline_envelope(Fi, tgrid);

    EbMat(s,:) = Eb;
    EiMat(s,:) = Ei;

    s_here = 100 * (Eb - Ei) ./ Eb;
    s_here(isnan(Eb) | isnan(Ei) | Eb <= 0) = NaN; % valid range only
    sav(s,:) = s_here;

    domRate(s) = dominance_rate(Fb, Fi) * 100;
end

meanSav = nanmean_local(sav, 1);
stdSav  = nanstd_local(sav, 0, 1);

Eb_mean = nanmean_local(EbMat, 1);
Ei_mean = nanmean_local(EiMat, 1);

% ---------------- HV & IGD (stored or estimated) ----------------
[HVb, IGDb, haveStoredB] = extract_HV_IGD(runs(isBase));
[HVi, IGDi, haveStoredI] = extract_HV_IGD(runs(isImp));
needEstimate = ~(haveStoredB && haveStoredI);

% BUG FIX #4: Note about HV normalization for cross-IS comparison
if haveStoredB && haveStoredI
    % HV/IGD stored during benchmark are normalized per-run against union PF
    % Therefore NOT directly comparable across different IS segments
    % (each segment has different PF reference)
end

if needEstimate
    % Estimate HV/IGD from fronts (normalized by union Pareto)
    Fall = collect_all_F(runs);
    Fref = nondominated_set(Fall);
    [HVb, IGDb] = estimate_hv_igd(runs(isBase), Fref);
    [HVi, IGDi] = estimate_hv_igd(runs(isImp),  Fref);
end

% ---------------- plots ----------------
plot_pareto(Fb_nd_all, Fi_nd_all, solver, nsga2_variant);
plot_saving_curve_union(tgrid, meanSav, stdSav, solver, nsga2_variant);     % Option 1
plot_deadline_envelopes(tgrid, Eb_mean, Ei_mean, solver, nsga2_variant);    % Option 2
plot_box_metric(HVb, HVi, 'HV - Hyper Volume  (higher is better)', solver, needEstimate, nsga2_variant);
plot_box_metric(IGDb, IGDi, 'IGD - Inverted Generational Distance (lower is better)', solver, needEstimate, nsga2_variant);

% ---------------- report ----------------
R = struct();
R.solver = solver;
R.tgrid  = tgrid;
R.meanSaving_percent = meanSav;
R.stdSaving_percent  = stdSav;
R.meanSaving_over_grid = nanmean_local(meanSav, 2);
R.minSaving_over_grid  = min(meanSav(~isnan(meanSav)));
R.maxSaving_over_grid  = max(meanSav(~isnan(meanSav)));
R.domRate_mean = nanmean_local(domRate, 1);
R.domRate_std  = nanstd_local(domRate, 0, 1);

% Store Pareto fronts for analysis scripts (needed for batch analysis)
R.F_b_union = Fb_nd_all;
R.F_i_union = Fi_nd_all;

fprintf('\n=== CCCR saving report | solver=%s ===\n', solver);
fprintf('Time axis shown (union): [%.2f, %.2f] s\n', tmin_union, tmax_union);
fprintf('Mean saving across shown axis: %.2f%%\n', R.meanSaving_over_grid);
fprintf('Min/Max mean saving across shown axis: %.2f%% / %.2f%%\n', R.minSaving_over_grid, R.maxSaving_over_grid);
fprintf('Dominance rate (base dominated by improved): %.1f ± %.1f%% (across seeds)\n', R.domRate_mean, R.domRate_std);
fprintf('HV: higher is better | IGD: lower is better\n');
if needEstimate
    fprintf('Note: HV/IGD were ESTIMATED (not stored in results).\n');
end

% ---------------- exports ----------------
T1 = table(tgrid(:), meanSav(:), stdSav(:), ...
    'VariableNames', {'T_deadline_s','Saving_mean_percent','Saving_std_percent'});
writetable(T1, sprintf('saving_curve_%s.xlsx', solver));

T2 = table(tgrid(:), Eb_mean(:), Ei_mean(:), ...
    'VariableNames', {'T_deadline_s','E_base_min_kWh','E_improved_min_kWh'});
writetable(T2, sprintf('deadline_energy_%s.xlsx', solver));

T3 = table([HVb(:); HVi(:)], [IGDb(:); IGDi(:)], ...
    [repmat("base",numel(HVb),1); repmat("improved",numel(HVi),1)], ...
    'VariableNames', {'HV','IGD','strategy'});
writetable(T3, sprintf('metrics_%s.xlsx', solver));

fprintf('Saved: saving_curve_%s.xlsx, deadline_energy_%s.xlsx, metrics_%s.xlsx\n', solver, solver, solver);

end

% ===================== helper functions =====================

function F = collect_all_F(runsSubset)
F = [];
for k = 1:numel(runsSubset)
    Fk = get_F_from_run(runsSubset(k));
    if ~isempty(Fk), F = [F; Fk]; end %#ok<AGROW>
end
end

function F = get_F_from_run(run)
% Try multiple common fields, then fallback to extract [T,E] from population.
cand = {'F','front','pareto','PF','paretoF','F_nd','F_pareto'};
F = [];

for i = 1:numel(cand)
    if isfield(run, cand{i})
        F = run.(cand{i});
        break;
    end
end

% fallback: if run contains "pop" (N x (D+4)) extract [T,E] correctly
if isempty(F)
    if isfield(run,'pop')
        F = run.pop;
    elseif isfield(run,'OUT') && isstruct(run.OUT) && isfield(run.OUT,'pop')
        F = run.OUT.pop;
    end
end

F = double(F);
if isempty(F) || ~ismatrix(F)
    F = [];
    return;
end

nc = size(F,2);
if nc == 2
    % already [T E] (or [E T])
elseif nc > 2
    % BUG FIX #1: Extract [T,E] CORRECTLY, not last 2 cols (which are rank/SDE)
    global dimension
    if ~isempty(dimension) && dimension > 0 && nc >= dimension + 2
        % Standard format: [X1...Xd | T | E | rank | SDE/crowd]
        F = F(:, dimension+1 : dimension+2);  % Extract T, E
    else
        % Fallback: assume [X... | T | E | rank | SDE]
        F = F(:, (nc-3):(nc-2));  % 2 cols before last 2
    end
else
    F = [];
    return;
end

% BUG FIX #2: Filter penalized solutions BEFORE column swap check
% (penalti solutions like T=1e6, E=1e6 cause axis scaling to zoom out)
MAX_T_REASONABLE = 1e5;  % threshold for penalized solutions
MAX_E_REASONABLE = 1e5;
valid = F(:,1) > 0 & F(:,1) < MAX_T_REASONABLE & ...
        F(:,2) >= 0 & F(:,2) < MAX_E_REASONABLE & ...
        all(isfinite(F), 2) & all(isreal(F), 2);
F = F(valid, :);
if isempty(F), return; end

% If columns are swapped (rare), try to detect and fix:
% Typical: Time ~ 50..400 s, Energy ~ 0..50 kWh
c1 = F(:,1); c2 = F(:,2);
q1 = quantile_safe(c1, 0.9);
q2 = quantile_safe(c2, 0.9);
if q1 < 60 && q2 > 60
    % likely [E T] -> swap to [T E]
    F = [c2 c1];
end
end

function q = quantile_safe(x, p)
x = x(~isnan(x));
if isempty(x), q = NaN; return; end
x = sort(x(:));
idx = max(1, min(numel(x), round(p*numel(x))));
q = x(idx);
end

function Fnd = nondominated_set(F)
if isempty(F), Fnd = zeros(0,2); return; end
keep = true(size(F,1),1);
for i=1:size(F,1)
    if ~keep(i), continue; end
    for j=1:size(F,1)
        if i==j || ~keep(j), continue; end
        if all(F(j,:)<=F(i,:)) && any(F(j,:)<F(i,:))
            keep(i)=false; break;
        end
    end
end
Fnd = F(keep,:);
Fnd = sortrows(Fnd,1,'ascend');
end

function E = deadline_envelope(F, tgrid)
% Min energy achievable under constraint T <= deadline
F = sortrows(F,1,'ascend');
cumMinE = cummin(F(:,2));
E = nan(size(tgrid));
k = 1;
for i=1:numel(tgrid)
    while k <= size(F,1) && F(k,1) <= tgrid(i) + 1e-9
        k = k + 1;
    end
    kk = k - 1;
    if kk >= 1
        E(i) = cumMinE(kk);
    end
end
end

function rate = dominance_rate(Fbase, Fimp)
if isempty(Fbase) || isempty(Fimp), rate = NaN; return; end
dom = false(size(Fbase,1),1);
for i=1:size(Fbase,1)
    b = Fbase(i,:);
    dom(i) = any(all(Fimp<=b,2) & any(Fimp<b,2));
end
rate = mean(dom);
end

function plot_pareto(Fb, Fi, solver, variant)
    if nargin < 4, variant = ''; end
    variantStr = '';
    if ~isempty(variant) && ~strcmp(variant, 'unknown')
        variantStr = sprintf(' [%s]', variant);
    end
    
    figure('Name', sprintf('Pareto (%s%s)', solver, variantStr));
    scatter(Fb(:,1), Fb(:,2), 28, [0.2 0.4 0.8], 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', 'base CC-CR');
    hold on;
    scatter(Fi(:,1), Fi(:,2), 20, [0.9 0.2 0.2], 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', 'improved CC-CR');
    grid on;
    xlabel('Time T (s)'); ylabel('Energy E (kWh)');
    title(sprintf('Pareto Front | %s%s | base=%d pts, imp=%d pts', ...
          solver, variantStr, size(Fb,1), size(Fi,1)));
    legend('Location','best');

    % Check if identical — size-safe comparison
    isIdentical = isequal(Fb, Fi);
    if ~isIdentical && isequal(size(Fb), size(Fi))
        isIdentical = max(max(abs(Fb - Fi))) < 1e-6;
    end
    if isIdentical
        text(0.05, 0.95, 'Note: base = improved (identical)', ...
             'Units','normalized', 'Color','red', 'FontSize',9);
    end
end

function plot_saving_curve_union(tgrid, meanSav, stdSav, solver, variant)
    if nargin < 5, variant = ''; end
    variantStr = '';
    if ~isempty(variant) && ~strcmp(variant, 'unknown')
        variantStr = sprintf(' [%s]', variant);
    end
    
    figure('Name', sprintf('Saving vs deadline (%s%s)', solver, variantStr));
    plot(tgrid, meanSav, 'LineWidth', 1.8); hold on;
    plot(tgrid, meanSav+stdSav, '--', 'LineWidth', 1.0);
    plot(tgrid, meanSav-stdSav, '--', 'LineWidth', 1.0);
    yline(0, ':', 'LineWidth', 1.2); % 0% line
    grid on;
    xlabel('Deadline T_{max} (s) (union axis)');
    ylabel('Saving (%) = (E_{base}-E_{imp})/E_{base} * 100');
    title(sprintf('Option 1: Energy Saving vs Deadline | %s%s', solver, variantStr));
    legend('mean','mean+std','mean-std','0% line','Location','best');
end

function plot_deadline_envelopes(tgrid, Eb_mean, Ei_mean, solver, variant)
    if nargin < 5, variant = ''; end
    variantStr = '';
    if ~isempty(variant) && ~strcmp(variant, 'unknown')
        variantStr = sprintf(' [%s]', variant);
    end
    
    figure('Name', sprintf('Deadline-energy envelopes (%s%s)', solver, variantStr));
    plot(tgrid, Eb_mean, 'LineWidth', 1.8); hold on;
    plot(tgrid, Ei_mean, 'LineWidth', 1.8);
    grid on;
    xlabel('Deadline T_{max} (s) (union axis)');
    ylabel('Min energy achievable (kWh) with T \le T_{max}');
    title(sprintf('Option 2: Deadline-energy envelopes | %s%s', solver, variantStr));
    legend('base(CC-CR)','improved (CC-CR)','Location','best');
end

function plot_box_metric(baseVals, impVals, ylab, solver, isEstimated, variant)
    if nargin < 6, variant = ''; end
    variantStr = '';
    if ~isempty(variant) && ~strcmp(variant, 'unknown')
        variantStr = sprintf(' [%s]', variant);
    end
    
    figure('Name', sprintf('%s (%s%s)', ylab, solver, variantStr));
    vals = [baseVals(:); impVals(:)];
    grp  = [repmat({'CC-CR'}, numel(baseVals), 1); repmat({'improved CC-CR'}, numel(impVals), 1)];
    boxplot(vals, grp);
    grid on;
    ttl = sprintf('%s | %s%s', ylab, solver, variantStr);
    if isEstimated
        ttl = [ttl ' (estimated)'];
    end
    title(ttl);
    ylabel(ylab);
end

function [HV, IGD, haveStored] = extract_HV_IGD(runsSubset)
HV = []; IGD = []; haveStored = true;

% Direct fields
if isfield(runsSubset,'HV'),  HV  = [runsSubset.HV];  end
if isfield(runsSubset,'IGD'), IGD = [runsSubset.IGD]; end

% Nested metrics struct
if (isempty(HV) || isempty(IGD)) && isfield(runsSubset,'metrics')
    try
        if isempty(HV),  HV  = arrayfun(@(r) r.metrics.HV,  runsSubset); end
        if isempty(IGD), IGD = arrayfun(@(r) r.metrics.IGD, runsSubset); end
    catch
        % ignore
    end
end

HV = HV(:); IGD = IGD(:);
haveStored = ~isempty(HV) && ~isempty(IGD);
end

function [HV, IGD] = estimate_hv_igd(runsSubset, Fref)
% Normalize using union Pareto range, then:
% - IGD: mean distance from each ref point to nearest point in run front
% - HV: 2D minimization hypervolume w.r.t ref point (1,1)
Fref = double(Fref);
lb = min(Fref,[],1);
ub = max(Fref,[],1);
rngv = max(ub-lb, 1e-9);
FrefN = (Fref - lb) ./ rngv;

HV = nan(numel(runsSubset),1);
IGD = nan(numel(runsSubset),1);

for k = 1:numel(runsSubset)
    Fk = get_F_from_run(runsSubset(k));
    Fk = nondominated_set(Fk);
    if isempty(Fk), continue; end
    FkN = (Fk - lb) ./ rngv;

    % IGD (use pdist2 if available)
    try
        d = pdist2(FrefN, FkN);
        IGD(k) = nanmean_local(min(d,[],2), 1);
    catch
        IGD(k) = igd_fallback(FrefN, FkN);
    end

    HV(k) = hypervolume2d(FkN, [1 1]);
end
end

function igd = igd_fallback(A, B)
% A: ref (m x 2), B: front (n x 2)
m = size(A,1);
mins = nan(m,1);
for i=1:m
    di = sqrt(sum((B - A(i,:)).^2, 2));
    mins(i) = min(di);
end
igd = nanmean_local(mins, 1);
end

function hv = hypervolume2d(F, ref)
% 2D HV for minimization
F = F(~any(isnan(F),2),:);
if isempty(F), hv = NaN; return; end
F = nondominated_set(F);

% clip
F(:,1) = min(F(:,1), ref(1));
F(:,2) = min(F(:,2), ref(2));

F = sortrows(F,1,'ascend');
hv = 0;
prevE = ref(2);
for i = 1:size(F,1)
    t = F(i,1);
    e = F(i,2);
    hv = hv + (ref(1) - t) * max(prevE - e, 0);
    prevE = min(prevE, e);
end
end

% ---- robust NaN-mean/std that won't break if "mean" is shadowed ----
function m = nanmean_local(A, dim)
if nargin < 2, dim = 1; end
mask = ~isnan(A);
cnt  = sum(mask, dim);

A(~mask) = 0;
m = sum(A, dim) ./ max(cnt, 1);

% IMPORTANT: if no data at all -> NaN (not 0)
m(cnt == 0) = NaN;
end

function s = nanstd_local(A, flag, dim)
if nargin < 2, flag = 0; end
if nargin < 3, dim = 1; end

mask = ~isnan(A);
cnt  = sum(mask, dim);

m = nanmean_local(A, dim);
X = bsxfun(@minus, A, m);
X(~mask) = 0;

den = cnt - (flag==0);   % unbiased if flag==0
den = max(den, 1);

s = sqrt(sum(X.^2, dim) ./ den);

% if no data -> NaN
s(cnt == 0) = NaN;
end
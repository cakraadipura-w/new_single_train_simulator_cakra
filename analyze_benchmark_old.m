function R = analyze_benchmark_old(matfile, solver, T_window, ngrid)
%ANALYZE_CCCR_BENCHMARK  Plot Pareto, HV/IGD, and % energy saving improved vs base.
%
% Usage:
%   R = analyze_cccr_benchmark('benchmark_results.mat','nsga2',[150 200],200);
%
% Notes:
% - Assumes objective points are F(:,1)=T [s], F(:,2)=E [kWh].
% - If HV/IGD are not stored in results, the script will estimate HV/IGD
%   using union Pareto as reference (normalized), so values may differ from yours.

if nargin < 1 || isempty(matfile), matfile = 'benchmark_results.mat'; end
if nargin < 2 || isempty(solver),  solver  = 'nsga2'; end
if nargin < 3 || isempty(T_window), T_window = [150 200]; end
if nargin < 4 || isempty(ngrid),   ngrid   = 200; end

if ~isfile(matfile)
    error('File not found: %s (check Current Folder / path)', matfile);
end

S = load(matfile,'results');
runs = S.results;

% ---- filter by solver ----
runs = runs(strcmp({runs.solver}, solver));
if isempty(runs), error('No runs found for solver="%s".', solver); end

isBase = strcmp({runs.strategy}, 'base');
isImp  = strcmp({runs.strategy}, 'improved');

if ~any(isBase) || ~any(isImp)
    error('Need both base and improved strategies in results.');
end

seedsBase = unique([runs(isBase).seed]);
seedsImp  = unique([runs(isImp).seed]);
seeds = intersect(seedsBase, seedsImp);
if isempty(seeds)
    error('No common seeds between base and improved.');
end

% ---- collect Pareto points ----
Fb_all = collect_all_F(runs(isBase));
Fi_all = collect_all_F(runs(isImp));
Fb_nd_all = nondominated_set(Fb_all);
Fi_nd_all = nondominated_set(Fi_all);

% overlap range
tmin = max(min(Fb_nd_all(:,1)), min(Fi_nd_all(:,1)));
tmax = min(max(Fb_nd_all(:,1)), max(Fi_nd_all(:,1)));

if ~isempty(T_window)
    tmin = max(tmin, T_window(1));
    tmax = min(tmax, T_window(2));
end
if tmax <= tmin
    error('No overlapping time range after applying T_window.');
end
tgrid = linspace(tmin, tmax, ngrid);

% ---- savings per seed ----
sav = nan(numel(seeds), ngrid);
Eb_mat = nan(numel(seeds), ngrid);
Ei_mat = nan(numel(seeds), ngrid);
domRate = nan(numel(seeds),1);

for s = 1:numel(seeds)
    sd = seeds(s);

    Fb = collect_all_F(runs(isBase & [runs.seed]==sd));
    Fi = collect_all_F(runs(isImp  & [runs.seed]==sd));
    Fb = nondominated_set(Fb);
    Fi = nondominated_set(Fi);

    Eb = deadline_envelope(Fb, tgrid); % min E with T<=deadline
    Ei = deadline_envelope(Fi, tgrid);

    Eb_mat(s,:) = Eb;
    Ei_mat(s,:) = Ei;

    sav(s,:) = 100 * (Eb - Ei) ./ Eb;
    domRate(s) = dominance_rate(Fb, Fi) * 100;
end

meanSav = nanmean_local(sav, 1);
stdSav  = nanstd_local(sav, 0, 1);

% ---- HV & IGD: use stored if available, else estimate ----
[HVb, IGDb, HVb_isStored, IGDb_isStored] = extract_metrics(runs(isBase));
[HVi, IGDi, HVi_isStored, IGDi_isStored] = extract_metrics(runs(isImp));

needEstimate = ~(HVb_isStored && HVi_isStored && IGDb_isStored && IGDi_isStored);
if needEstimate
    % Estimate HV/IGD using union Pareto reference and normalization
    Fall = collect_all_F(runs);
    Fref = nondominated_set(Fall);

    [HVb, IGDb] = estimate_hv_igd(runs(isBase), Fref);
    [HVi, IGDi] = estimate_hv_igd(runs(isImp),  Fref);
end

% ---- plots ----
plot_pareto(Fb_nd_all, Fi_nd_all, solver);
plot_saving_curve(tgrid, meanSav, stdSav, solver);
plot_box_metrics(HVb, HVi, 'HV', solver, needEstimate);
plot_box_metrics(IGDb, IGDi, 'IGD', solver, needEstimate);

% ---- report + export ----
R = struct();
R.solver = solver;
R.tgrid = tgrid;
R.meanSaving_percent = meanSav;
R.stdSaving_percent  = stdSav;
R.meanSaving_over_window = mean(meanSav,'omitnan');
R.medianSaving_over_window = median(meanSav,'omitnan');
R.minSaving_over_window = min(meanSav);
R.maxSaving_over_window = max(meanSav);
R.domRate_percent_mean = mean(domRate,'omitnan');
R.domRate_percent_std  = std(domRate,'omitnan');

fprintf('\n=== CCCR saving report | solver=%s ===\n', solver);
fprintf('Time overlap used: [%.2f, %.2f] s\n', tmin, tmax);
fprintf('Mean saving (improved vs base): %.2f%%  (median %.2f%%)\n', R.meanSaving_over_window, R.medianSaving_over_window);
fprintf('Min/Max mean saving in window: %.2f%% / %.2f%%\n', R.minSaving_over_window, R.maxSaving_over_window);
fprintf('Dominance rate (base dominated by improved): %.1f ± %.1f%% (across seeds)\n', R.domRate_percent_mean, R.domRate_percent_std);

% export saving curve
T1 = table(tgrid(:), meanSav(:), stdSav(:), 'VariableNames', {'T_deadline_s','Saving_mean_percent','Saving_std_percent'});
%writetable(T1, sprintf('saving_curve_%s.xlsx', solver));

% export per-seed envelope (optional, big but useful)
Eb_mean = nanmean_local(Eb_mat, 1);
Ei_mean = nanmean_local(Ei_mat, 1);

T2 = table(tgrid(:), Eb_mean(:), Ei_mean(:), 'VariableNames', {'T_deadline_s','E_base_min_kWh','E_improved_min_kWh'});
%writetable(T2, sprintf('deadline_energy_%s.xlsx', solver));

% export metrics
T3 = table(HVb(:), IGDb(:), repmat("base",numel(HVb),1), 'VariableNames', {'HV','IGD','strategy'});
T4 = table(HVi(:), IGDi(:), repmat("improved",numel(HVi),1), 'VariableNames', {'HV','IGD','strategy'});
%writetable([T3;T4], sprintf('metrics_%s.xlsx', solver));

fprintf('Saved: saving_curve_%s.xlsx, deadline_energy_%s.xlsx, metrics_%s.xlsx\n', solver, solver, solver);

end

% ===================== helpers =====================

function F = collect_all_F(runsSubset)
F = [];
for k = 1:numel(runsSubset)
    Fk = get_F_from_run(runsSubset(k));
    if ~isempty(Fk), F = [F; Fk]; end %#ok<AGROW>
end
end

function F = get_F_from_run(run)
% Try common field names
cand = {'F','front','pareto','PF','paretoF','F_nd','F_pareto'};
F = [];
for i = 1:numel(cand)
    if isfield(run, cand{i})
        F = run.(cand{i});
        break;
    end
end
if isempty(F)
    % Try nested
    if isfield(run,'OUT') && isstruct(run.OUT) && isfield(run.OUT,'F')
        F = run.OUT.F;
    elseif isfield(run,'metrics') && isstruct(run.metrics) && isfield(run.metrics,'F')
        F = run.metrics.F;
    end
end
if isempty(F), return; end

F = double(F);

if isempty(F) || ~ismatrix(F)
    F = [];
    return;
end

nc = size(F,2);
if nc == 2
    % OK: already [T E]
elseif nc > 2
    F = F(:, (nc-1):nc);     % <-- aman, nggak pakai end-1
else
    F = [];                  % nc < 2, nggak bisa bikin [T E]
    return;
end

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
% best energy under deadline constraint T <= t
F = sortrows(F,1,'ascend');
cumMinE = cummin(F(:,2));
E = nan(size(tgrid));
k = 1;
for i=1:numel(tgrid)
    while k <= size(F,1) && F(k,1) <= tgrid(i) + 1e-9
        k = k + 1;
    end
    kk = k-1;
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

function plot_pareto(Fb, Fi, solver)
figure('Name', sprintf('Pareto (%s)', solver));
scatter(Fb(:,1), Fb(:,2), 22, 'filled'); hold on;
scatter(Fi(:,1), Fi(:,2), 22, 'filled');
grid on;
xlabel('Time T (s)'); ylabel('Energy E (kWh)');
title(sprintf('Pareto Front Comparison | %s', solver));
legend('base','improved','Location','best');
end

function plot_saving_curve(tgrid, meanSav, stdSav, solver)
figure('Name', sprintf('Saving curve (%s)', solver));
plot(tgrid, meanSav, 'LineWidth', 1.8); hold on;
plot(tgrid, meanSav+stdSav, '--', 'LineWidth', 1.0);
plot(tgrid, meanSav-stdSav, '--', 'LineWidth', 1.0);
grid on;
xlabel('Deadline T_{max} (s)');
ylabel('Energy saving (%) = (E_{base}-E_{imp})/E_{base} * 100');
title(sprintf('Energy Saving vs Deadline | %s', solver));
legend('mean','mean+std','mean-std','Location','best');
end

function plot_box_metrics(HVb, HVi, metricName, solver, isEstimated)
figure('Name', sprintf('%s (%s)', metricName, solver));
vals = [HVb(:); HVi(:)];
grp  = [repmat({'base'}, numel(HVb), 1); repmat({'improved'}, numel(HVi), 1)];
boxplot(vals, grp);
grid on;
ttl = sprintf('%s Comparison | %s', metricName, solver);
if isEstimated
    ttl = [ttl ' (estimated)'];
end
title(ttl);
ylabel(metricName);
end

function [HV, IGD, hvStored, igdStored] = extract_metrics(runsSubset)
HV = []; IGD = [];
hvStored = true; igdStored = true;

% Try direct fields
if isfield(runsSubset,'HV'),  HV = [runsSubset.HV];  end
if isfield(runsSubset,'IGD'), IGD = [runsSubset.IGD]; end

% Try nested metrics
if isempty(HV)
    hvStored = false;
    if isfield(runsSubset,'metrics')
        try
            HV = arrayfun(@(r) r.metrics.HV, runsSubset);
            hvStored = true;
        catch
            HV = [];
        end
    end
end
if isempty(IGD)
    igdStored = false;
    if isfield(runsSubset,'metrics')
        try
            IGD = arrayfun(@(r) r.metrics.IGD, runsSubset);
            igdStored = true;
        catch
            IGD = [];
        end
    end
end

if isempty(HV),  hvStored = false; end
if isempty(IGD), igdStored = false; end

HV = HV(:);
IGD = IGD(:);
end

function [HV, IGD] = estimate_hv_igd(runsSubset, Fref)
% Normalize by union range and compute:
% - IGD: mean dist from each ref point to nearest point in run front
% - HV: 2D minimization hypervolume (normalized) using ref point (1,1)

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

    % IGD
    d = pdist2(FrefN, FkN);
    IGD(k) = mean(min(d,[],2), 'omitnan');

    % HV in 2D (minimization), reference point at (1,1)
    HV(k) = hypervolume2d(FkN, [1 1]);
end
end

function hv = hypervolume2d(F, ref)
% Compute 2D hypervolume for minimization.
% Assumes F is nondominated in minimization sense.
F = F(~any(isnan(F),2),:);
if isempty(F), hv = NaN; return; end
F = nondominated_set(F);
% Clip to ref
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

function m = nanmean_local(A, dim)
if nargin < 2, dim = 1; end
mask = ~isnan(A);
A(~mask) = 0;
m = sum(A, dim) ./ max(1, sum(mask, dim));
end

function s = nanstd_local(A, flag, dim)
if nargin < 2, flag = 0; end
if nargin < 3, dim = 1; end

m = nanmean_local(A, dim);
X = bsxfun(@minus, A, m);

mask = ~isnan(A);
X(~mask) = 0;

n = sum(mask, dim);
den = n - (flag==0);        % unbiased if flag==0
den = max(den, 1);

s = sqrt(sum(X.^2, dim) ./ den);
end
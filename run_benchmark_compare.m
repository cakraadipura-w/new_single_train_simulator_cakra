function [results, compare] = run_benchmark_compare(cfg)
%RUN_BENCHMARK_COMPARE  Benchmark solvers across seeds AND driving strategies.
%
% Outputs:
%   results : struct array (one entry per run: strategy x solver x seed)
%             Includes AllX/AllF (final population) + Pareto X/F for that run.
%   compare : struct array (one entry per solver) with energy-saving (%) tables
%             comparing improved vs base along time limits.
%
% Assumptions:
%   - Objectives are MINIMIZED: F = [Time_s, Energy_kWh]
%   - Solver returns a numeric population matrix, typically:
%       pop = [X(1:dim), f1, f2, rank, crowd]  (NSGA-II style)
%
% Required globals (usually set by setup_project / init_worker_globals):
%   global pop_size iterations dimension vel_profile use_improved
%
% If cfg.parallel_use=true and a pool exists, this function calls:
%   parfevalOnAll(@init_worker_globals, 0, routePath, rollingstockPath, cfg.sim, use_improved)
% whenever strategy changes, so workers refresh globals.

% ------------------- Globals -------------------
global pop_size iterations dimension vel_profile
global use_improved

% ------------------- Config -------------------
solvers    = lower(string(cfg.bench_solvers));
seeds      = cfg.bench_seeds;
make_plots = isfield(cfg,'make_plots') && cfg.make_plots;

% Strategies: {'base','improved'} recommended
if isfield(cfg,'bench_strategies') && ~isempty(cfg.bench_strategies)
    strategies = lower(string(cfg.bench_strategies));
else
    strategies = string(use_improved);
end

% Paths for worker refresh
if isfield(cfg,'routePath') && ~isempty(cfg.routePath)
    routePath = cfg.routePath;
else
    routePath = which(cfg.route_file);
end
if isfield(cfg,'rollingstockPath') && ~isempty(cfg.rollingstockPath)
    rollingstockPath = cfg.rollingstockPath;
else
    rollingstockPath = which(cfg.rollingstock_file);
end

% Fair iteration budget (optional override)
benchPop = cfg.pop_size;
benchIters = struct();
for s = 1:numel(solvers)
    benchIters.(char(solvers(s))) = cfg.iterations;
end
if isfield(cfg,'bench_iters') && ~isempty(cfg.bench_iters)
    f = fieldnames(cfg.bench_iters);
    for i = 1:numel(f)
        benchIters.(lower(f{i})) = cfg.bench_iters.(f{i});
    end
end

% Backup globals so we can restore
pop0 = pop_size;
it0  = iterations;

% Template (fixed fields)
tmpl = struct( ...
    'strategy','', ...
    'solver','', ...
    'seed',0, ...
    'pop_size',0, ...
    'dim',0, ...
    'iterations',0, ...
    'runtime',NaN, ...
    'AllX',zeros(0,0), ...
    'AllF',zeros(0,2), ...
    'X',zeros(0,0), ...
    'F',zeros(0,2), ...
    'HV',NaN, ...
    'IGD',NaN, ...
    'Nf1',0);

runs = repmat(tmpl, 0, 1);
allFronts = zeros(0,2);

% ------------------- Run benchmark -------------------
fprintf('\nRunning benchmark: strategies=%s | solvers=%s | seeds=%s\n', ...
    strjoin(cellstr(strategies),','), strjoin(cellstr(solvers),','), mat2str(seeds));

for st = 1:numel(strategies)
    strat = char(strategies(st));

    % Switch global use_improved based on strategy name
    if strcmpi(strat,'improved') || strcmpi(strat,'1') || strcmpi(strat,'true')
        use_improved = true;
    else
        use_improved = false;
    end

    % Refresh globals on client (ensures dimension/vel_profile consistent per strategy)
    if exist('init_worker_globals','file') == 2
        if ~isfield(cfg,'sim'), cfg.sim = struct(); end
        init_worker_globals(routePath, rollingstockPath, cfg.sim, use_improved);
    end

    fprintf('\n--- Strategy=%s | use_improved=%d ---\n', strat, use_improved);

    % Refresh worker globals when using parallel pool
    if isfield(cfg,'parallel_use') && cfg.parallel_use
        pool = gcp('nocreate');
        if ~isempty(pool) && exist('init_worker_globals','file') == 2
            if ~isfield(cfg,'sim'), cfg.sim = struct(); end
            fInit = parfevalOnAll(@init_worker_globals, 0, routePath, rollingstockPath, cfg.sim, use_improved);
            fetchOutputs(fInit);
        end
    end

    for s = 1:numel(solvers)
        solver = solvers(s);

        pop_size   = benchPop;
        iterations = benchIters.(char(solver));

        for k = 1:numel(seeds)
            rng(seeds(k),'twister');

            fprintf('Running %s | %s | seed=%d | pop=%d | iter=%d\n', strat, solver, seeds(k), pop_size, iterations);

            t0 = tic;
            pop = run_one_solver(solver, vel_profile);
            rt = toc(t0);

            % Extract ALL population points (final generation)
            [AllX, AllF] = extract_all_XF(pop, dimension);

            % Extract Pareto (nondominated) + matching decision variables
            [PX, PF] = extract_pareto_XF(pop, dimension, AllX, AllF);

            r = tmpl;
            r.strategy   = strat;
            r.solver     = char(solver);
            r.seed       = seeds(k);
            r.pop_size   = size(AllF,1);
            r.dim        = size(AllX,2);
            r.iterations = iterations;
            r.runtime    = rt;

            r.AllX = AllX;
            r.AllF = AllF;
            r.X    = PX;
            r.F    = PF;
            r.Nf1  = size(PF,1);

            runs(end+1,1) = r; %#ok<AGROW>

            if ~isempty(PF)
                allFronts = [allFronts; PF]; %#ok<AGROW>
            end
        end
    end
end

% Restore
pop_size  = pop0;
iterations = it0;

% ------------------- Metrics (HV/IGD) -------------------
compare = struct([]); % default output if early return

if isempty(allFronts)
    warning('No fronts collected. Check solver output / objective extraction.');
    results = runs;
    return;
end

% Build reference PF from union
PF_ref = nondominated_set(allFronts);

% Normalization bounds (robust: ignore crazy penalties)
mins = min(PF_ref, [], 1);
maxs = max(PF_ref, [], 1);
ref  = maxs + 0.10*(maxs - mins + eps); % 10% worse

Rn = (PF_ref - mins) ./ (maxs - mins + eps);
refn = (ref - mins) ./ (maxs - mins + eps);

for i = 1:numel(runs)
    F = runs(i).F;
    if isempty(F)
        runs(i).HV = 0;
        runs(i).IGD = NaN;
        runs(i).Nf1 = 0;
        continue;
    end
    Fn = (F - mins) ./ (maxs - mins + eps);
    runs(i).HV  = hv2d(Fn, refn);
    runs(i).IGD = igd(Rn, Fn);
end

% ------------------- Summary -------------------
fprintf('\n==== SUMMARY (mean ± std) ====\n');
for st = 1:numel(strategies)
    strat = char(strategies(st));
    for s = 1:numel(solvers)
        sol = char(solvers(s));
        idx = strcmpi({runs.strategy}, strat) & strcmpi({runs.solver}, sol);
        if ~any(idx), continue; end

        HV  = [runs(idx).HV];
        IGD = [runs(idx).IGD];
        RT  = [runs(idx).runtime];
        Nf1 = [runs(idx).Nf1];

        fprintf('%-8s | %-6s | HV: %.4f±%.4f | IGD: %.4f±%.4f | RT: %.1f±%.1f s | F1: %.1f±%.1f\n', ...
            strat, sol, mean(HV,'omitnan'), std(HV,'omitnan'), ...
            mean(IGD,'omitnan'), std(IGD,'omitnan'), ...
            mean(RT,'omitnan'), std(RT,'omitnan'), ...
            mean(Nf1,'omitnan'), std(Nf1,'omitnan'));
    end
end

% ------------------- Energy saving (%) compare -------------------
compare = compute_energy_compare(runs, solvers, strategies);

% Save convenience outputs (optional; comment out if you don't want files)
try
    save('energy_compare.mat','compare');
    for ii = 1:numel(compare)
        if isfield(compare(ii),'grid') && ~isempty(compare(ii).grid)
            writetable(compare(ii).grid, sprintf('energy_compare_%s.csv', compare(ii).solver));
        end
        if isfield(compare(ii),'at_basePF') && ~isempty(compare(ii).at_basePF)
            writetable(compare(ii).at_basePF, sprintf('energy_compare_at_basePF_%s.csv', compare(ii).solver));
        end
        if isfield(compare(ii),'at_improvedPF') && ~isempty(compare(ii).at_improvedPF)
            writetable(compare(ii).at_improvedPF, sprintf('energy_compare_at_improvedPF_%s.csv', compare(ii).solver));
        end
    end
catch ME
    warning("Failed to save energy_compare outputs: %s", ME.message);
end

% ------------------- Plots -------------------
if make_plots
    plot_pareto_overlay(runs, solvers, strategies);
    plot_energy_compare(compare);
end

results = runs;
end

%% ======================= helpers =======================

function pop = run_one_solver(solver, vel_profile)
% Wrapper to call your solvers (add more as needed)
switch char(solver)
    case 'nsga2'
        [pop,~,~] = nsga2_main(vel_profile);
    case 'mopso'
        if exist('mopso_main','file')~=2
            error('mopso_main.m not found.');
        end
        pop = mopso_main(vel_profile);
    case 'spea2'
        if exist('spea2_main','file')~=2
            error('spea2_main.m not found.');
        end
        pop = spea2_main(vel_profile);
    case 'moead'
        if exist('moead_main','file')~=2
            error('moead_main.m not found.');
        end
        pop = moead_main(vel_profile);
    otherwise
        error('Unknown solver: %s', solver);
end
end

function [AllX, AllF] = extract_all_XF(pop, dim)
% Extract full population decision vars and objectives from solver output.
if isempty(pop)
    AllX = zeros(0,dim);
    AllF = zeros(0,2);
    return;
end
if ~isnumeric(pop)
    error('Expected numeric pop matrix. Got: %s', class(pop));
end
if size(pop,2) < dim+2
    error('pop has too few columns to contain X + [T E]. size=%dx%d dim=%d', size(pop,1), size(pop,2), dim);
end
AllX = double(pop(:, 1:dim));
AllF = double(pop(:, dim+1:dim+2)); % [T E]
end

function [PX, PF] = extract_pareto_XF(pop, dim, AllX, AllF)
% Return nondominated set + matching decision vectors.
if isempty(AllF)
    PX = zeros(0,dim);
    PF = zeros(0,2);
    return;
end

% Candidate points: use rank==1 if available (NSGA-II)
cand = (1:size(AllF,1))';
if isnumeric(pop) && size(pop,2) >= dim+4
    rankcol = dim+3;
    r = pop(:, rankcol);
    c = find(r==1);
    if ~isempty(c), cand = c; end
end

mask = nondominated_mask(AllF(cand,:));
idx  = cand(mask);

PF = AllF(idx,:);
PX = AllX(idx,:);

[~,ord] = sort(PF(:,1),'ascend');
PF = PF(ord,:);
PX = PX(ord,:);
end

function keep = nondominated_mask(F)
% 2D minimization nondominated logical mask
n = size(F,1);
keep = true(n,1);
for i = 1:n
    if ~keep(i), continue; end
    for j = 1:n
        if i==j || ~keep(j), continue; end
        if all(F(j,:)<=F(i,:)) && any(F(j,:)<F(i,:))
            keep(i) = false;
            break;
        end
    end
end
end

function Fnd = nondominated_set(F)
% Return nondominated set (Nx2) sorted by Time ascending
if isempty(F)
    Fnd = zeros(0,2);
    return;
end
keep = nondominated_mask(F);
Fnd = F(keep,:);
[~,idx] = sort(Fnd(:,1),'ascend');
Fnd = Fnd(idx,:);
end

function hv = hv2d(F, ref)
% 2D hypervolume for minimization in normalized space.
% ref must be worse than all points.
if isempty(F)
    hv = 0;
    return;
end
F = nondominated_set(F);

% enforce monotonic in E for stable area
E = cummin(F(:,2));
F(:,2) = E;

hv = 0;
prevE = ref(2);
for i = 1:size(F,1)
    t = F(i,1); e = F(i,2);
    if e > prevE, continue; end
    hv = hv + (ref(1) - t) * (prevE - e);
    prevE = e;
end
end

function val = igd(R, A)
% Inverted Generational Distance: avg distance from each ref point to approximation set
if isempty(R) || isempty(A)
    val = NaN;
    return;
end
d = zeros(size(R,1),1);
for i = 1:size(R,1)
    dr = A - R(i,:);
    d(i) = min(sqrt(sum(dr.^2,2)));
end
val = mean(d);
end

function compare = compute_energy_compare(runs, solvers, strategies)
% Build energy-saving (%) curves per solver comparing improved vs base.
strategies = lower(string(strategies));
baseName = "base"; impName = "improved";
if ~any(strategies==baseName) || ~any(strategies==impName)
    if numel(strategies) >= 2
        baseName = strategies(1);
        impName  = strategies(2);
    else
        baseName = strategies(1);
        impName  = strategies(1);
    end
end

compare = repmat(struct( ...
    'solver','', ...
    'base',char(baseName), ...
    'improved',char(impName), ...
    'env_base',struct('T',[],'Ebest',[]), ...
    'env_improved',struct('T',[],'Ebest',[]), ...
    'grid',table(), ...
    'at_basePF',table(), ...
    'at_improvedPF',table(), ...
    'summary',struct() ), numel(solvers), 1);

for si = 1:numel(solvers)
    sol = char(solvers(si));

    Fb = collect_union_pf(runs, baseName, sol);
    Fi = collect_union_pf(runs, impName,  sol);

    envB = make_energy_envelope(Fb);
    envI = make_energy_envelope(Fi);

    compare(si).solver = sol;
    compare(si).env_base = envB;
    compare(si).env_improved = envI;

    if isempty(envB.T) || isempty(envI.T)
        compare(si).summary = struct('note',"Missing data for one strategy");
        continue;
    end

    Tmin = max([min(envB.T), min(envI.T)]);
    Tmax = min([max(envB.T), max(envI.T)]);
    if Tmin < Tmax
        Tgrid = linspace(Tmin, Tmax, 80)';
    else
        Tgrid = unique([envB.T; envI.T]);
    end

    Eb = interp1(envB.T, envB.Ebest, Tgrid, 'pchip', 'extrap');
    Ei = interp1(envI.T, envI.Ebest, Tgrid, 'pchip', 'extrap');
    saving_pct = 100*(Eb - Ei)./max(Eb, eps);

    compare(si).grid = table(Tgrid, Eb, Ei, saving_pct, ...
        'VariableNames', {'T_limit_s','E_base_kWh','E_improved_kWh','saving_pct'});

    % "Per-Pareto-point" saving tables
    Tb = envB.T;
    Eb2 = envB.Ebest;
    Ei2 = interp1(envI.T, envI.Ebest, Tb, 'pchip', 'extrap');
    saveTb = 100*(Eb2 - Ei2)./max(Eb2, eps);
    compare(si).at_basePF = table(Tb, Eb2, Ei2, saveTb, ...
        'VariableNames', {'T_limit_s','E_base_kWh','E_improved_kWh','saving_pct'});

    Ti = envI.T;
    Ei3 = envI.Ebest;
    Eb3 = interp1(envB.T, envB.Ebest, Ti, 'pchip', 'extrap');
    saveTi = 100*(Eb3 - Ei3)./max(Eb3, eps);
    compare(si).at_improvedPF = table(Ti, Eb3, Ei3, saveTi, ...
        'VariableNames', {'T_limit_s','E_base_kWh','E_improved_kWh','saving_pct'});

    compare(si).summary = struct( ...
        'overlap_Tmin', Tmin, ...
        'overlap_Tmax', Tmax, ...
        'mean_saving_pct', mean(saving_pct,'omitnan'), ...
        'median_saving_pct', median(saving_pct,'omitnan'), ...
        'max_saving_pct', max(saving_pct), ...
        'min_saving_pct', min(saving_pct) );
end
end

function F = collect_union_pf(runs, strategyName, solverName)
% Collect union PF points across seeds for given strategy+solver
idx = strcmpi({runs.strategy}, char(strategyName)) & strcmpi({runs.solver}, char(solverName));
if ~any(idx)
    F = zeros(0,2);
    return;
end
F = vertcat(runs(idx).F);
F = double(F);
F = F(all(isfinite(F),2) & F(:,1)>0 & F(:,2)>=0, :);
F = unique(round(F,8), 'rows');
F = F(nondominated_mask(F),:);
F = sortrows(F,1);
end

function env = make_energy_envelope(F)
% Convert PF points [T,E] into best-energy-under-time-limit curve E*(Tlim)
if isempty(F)
    env = struct('T',zeros(0,1),'Ebest',zeros(0,1));
    return;
end
T = F(:,1); E = F(:,2);
[Ts, ord] = sort(T, 'ascend');
Es = E(ord);

Ebest = zeros(size(Es));
m = inf;
for i = 1:numel(Es)
    if Es(i) < m, m = Es(i); end
    Ebest(i) = m;
end

[Tu, ia] = unique(Ts,'stable');
Eu = Ebest(ia);
env = struct('T',Tu,'Ebest',Eu);
end

function plot_pareto_overlay(runs, solvers, strategies)
% One figure per solver: overlay union PF for each strategy
strategies = lower(string(strategies));
for si = 1:numel(solvers)
    sol = char(solvers(si));
    figure('Name', sprintf('Pareto overlay - %s', sol)); hold on; grid on;

    for st = 1:numel(strategies)
        strat = char(strategies(st));
        F = collect_union_pf(runs, strat, sol);
        if isempty(F), continue; end
        plot(F(:,1), F(:,2), '-o', 'LineWidth', 1.5, 'DisplayName', strat);
    end

    xlabel('Time T (s)'); ylabel('Energy E (kWh)');
    title(sprintf('Union Pareto Front overlay | solver=%s', sol));
    legend('Location','best');
end
end

function plot_energy_compare(compare)
% One solver per pair of figures: envelope + saving%
for i = 1:numel(compare)
    if isempty(compare(i).grid) || isempty(compare(i).env_base.T) || isempty(compare(i).env_improved.T)
        continue;
    end

    figure('Name', sprintf('E*(Tlim) - %s', compare(i).solver)); grid on; hold on;
    plot(compare(i).env_base.T, compare(i).env_base.Ebest, '-o', 'DisplayName', compare(i).base);
    plot(compare(i).env_improved.T, compare(i).env_improved.Ebest, '-o', 'DisplayName', compare(i).improved);
    xlabel('Time limit T_{lim} (s)'); ylabel('Best energy E*(Tlim) (kWh)');
    title(sprintf('Best energy under time limit | solver=%s', compare(i).solver));
    legend('Location','best');

    figure('Name', sprintf('Energy saving %% - %s', compare(i).solver)); grid on; hold on;
    plot(compare(i).grid.T_limit_s, compare(i).grid.saving_pct, 'LineWidth', 1.5);
    xlabel('Time limit T_{lim} (s)');
    ylabel('Energy saving (%)  (positive = improved saves)');
    title(sprintf('Saving = 100*(E_base - E_improved)/E_base | solver=%s', compare(i).solver));
end
end

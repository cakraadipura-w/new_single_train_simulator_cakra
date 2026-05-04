%% run_E1_ablation.m — Experiment E1: Ablation Study (G2)
%
% Membandingkan 4 kombinasi driving strategy × solver:
%   A — original CC_CR  + NSGA-II vanilla        (Baseline)
%   B — original CC_CR  + BRL-SDE-NSGA-II        (Solver only)
%   C — improved CC_CR  + NSGA-II vanilla        (Strategy only)
%   D — improved CC_CR  + BRL-SDE-NSGA-II        (Both, proposed)
%
% Pendekatan: seperti run_benchmark_compare.m
%   1. Jalankan setiap config N_RUNS kali
%   2. Bangun union Pareto per config (gabung semua runs)
%   3. Satu global reference front → HV/IGD comparable antar config
%   4. Grafik overlay Pareto + kurva energy saving (seperti plot_energy_compare)
%   5. Wilcoxon signrank test
%
% Outputs: E1_results.csv, E1_summary.txt,
%          E1_pareto_overlay.png, E1_energy_comparison.png,
%          E1_energy_combined.png, E1_hv_bar.png

clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% ===== RUN SWITCHES =====
% Set each flag to true to run that config, false to load from cache.
% Set RUN_ALL = true to override all switches and run every config.
%
%   Example — only re-run Config D:
%       RUN_CONFIG_A = false;
%       RUN_CONFIG_B = false;
%       RUN_CONFIG_C = false;
%       RUN_CONFIG_D = true;
%       RUN_ALL      = false;

RUN_CONFIG_A = true;
RUN_CONFIG_B = true;
RUN_CONFIG_C = true;
RUN_CONFIG_D = true;
RUN_ALL      = false;   % true = force-run ALL regardless of individual switches

%% ===== PER-CONFIG BUDGET =====
% All 0 = use global defaults (POP_SIZE / ITERATIONS / N_RUNS).
% Same budget for all configs — fair HV/IGD comparison (Option 2: alpha fixed).
%
%   CONFIG_POP  — population size per config  (0 = use POP_SIZE)
%   CONFIG_GEN  — number of generations       (0 = use ITERATIONS)
%   CONFIG_RUNS — number of independent runs  (0 = use N_RUNS)

CONFIG_POP  = struct('A', 500,   'B', 500,   'C', 500,   'D', 500);
CONFIG_GEN  = struct('A', 600,   'B', 660,   'C', 600,   'D', 600);
CONFIG_RUNS = struct('A', 1,     'B', 1,     'C', 1,     'D', 1);

%% ===== GLOBAL DEFAULTS =====
if exist('ACTIVE_SEG','var') && ~isempty(ACTIVE_SEG)
    ROUTE_IS04 = ACTIVE_SEG.file;
    T_TARGET   = ACTIVE_SEG.T_sched;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', ACTIVE_SEG.name);
    N_RUNS     = ACTIVE_NRUNS;
    SEG_LABEL  = ACTIVE_SEG.name;
else
    ROUTE_IS04 = 'Guangzhou_Line7_IS02_3.548-4.908km.mat';
    T_TARGET   = 170;
    N_RUNS     = 30;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', 'IS02');
    SEG_LABEL  = 'IS02';
end
RS_FILE    = 'rollingstock_Guangzhou_L7.m';
POP_SIZE   = 300;
ITERATIONS = 500;
N_WORKERS  = 4;

% Search bound: wider than T_TARGET so relaxed-schedule solutions appear in Pareto
SEARCH_TLIM = T_TARGET * 1.50;

% x-axis limit for Pareto plots
PLOT_TLIM   = T_TARGET * 1.30;

CONFIGS = struct( ...
    'id',            {'A',         'B',          'C',              'D'}, ...
    'use_improved',  {false,       false,         true,             true}, ...
    'nsga2_variant', {'original',  'rl_sde',     'original',       'rl_sde'}, ...
    'label',         {'Baseline',  'Solver only', 'Strategy only',  'Both (proposed)'}, ...
    'desc',          {'Original CC_CR + Original NSGA-II', ...
                      'Original CC_CR + BRL-SDE-NSGA-II', ...
                      'Improved CC_CR + Original NSGA-II', ...
                      'Improved CC_CR + BRL-SDE-NSGA-II'} );

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ===== GLOBALS + SETUP =====
global use_improved nsga2_variant pop_size iterations dimension
global vel_profile time_obj_max parallel_use show_progress driving_strategy

show_progress    = false;
parallel_use     = true;
driving_strategy = "CC_CR";
time_obj_max     = SEARCH_TLIM;   % lebih longgar — mode santai tetap bisa ditemukan

cfg = struct('route_file', ROUTE_IS04, 'rollingstock_file', RS_FILE, ...
    'driving_strategy', "CC_CR", 'pop_size', POP_SIZE, 'iterations', ITERATIONS, ...
    'use_improved', false, 'nsga2_variant', 'original', ...
    'parallel_use', true, 'sim', struct());
info = setup_project(cfg);

[~, E_BASELINE] = flat_out_baseline_rk4(info.routePath, info.rollingstockPath);
fprintf('E_baseline=%.4f kWh | T_target=%d s | Segment=%s\n', E_BASELINE, T_TARGET, SEG_LABEL);

setup_parallel_pool(cfg, info);
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', N_WORKERS);
end

%% ===== PHASE 1: RUN OR LOAD EACH CONFIG =====
n_cfg    = numel(CONFIGS);
% Size pops_all to the maximum runs any config might use
max_runs = N_RUNS;
for ci = 1:n_cfg
    r = CONFIG_RUNS.(CONFIGS(ci).id);
    if r > 0 && r > max_runs, max_runs = r; end
end
pops_all     = cell(n_cfg, max_runs);
t_mat        = nan(n_cfg, max_runs);
dim_all      = zeros(n_cfg, 1);
runs_per_cfg = zeros(n_cfg, 1);
pop_cfg      = zeros(n_cfg, 1);   % actual pop size used per config
gen_cfg      = zeros(n_cfg, 1);   % actual generations used per config

% Print run plan
fprintf('\n===== E1 Run Plan =====\n');
mat_cache = fullfile(OUT_DIR, 'E1_all_configs.mat');
for ci = 1:n_cfg
    id     = CONFIGS(ci).id;
    do_run = RUN_ALL || eval(['RUN_CONFIG_' id]);
    pop_ci  = CONFIG_POP.(id);  if pop_ci  == 0, pop_ci  = POP_SIZE;   end
    gen_ci  = CONFIG_GEN.(id);  if gen_ci  == 0, gen_ci  = ITERATIONS; end
    runs_ci = CONFIG_RUNS.(id); if runs_ci == 0, runs_ci = N_RUNS;     end
    if do_run
        fprintf('  Config %s — RUN  | pop=%d | gen=%d | runs=%d\n', id, pop_ci, gen_ci, runs_ci);
    else
        cached = exist(mat_cache,'file') == 2;
        fprintf('  Config %s — SKIP | %s\n', id, ternary(cached,'will load from cache','NO CACHE — results will be empty!'));
    end
end
fprintf('=======================\n\n');

for ci = 1:n_cfg
    id     = CONFIGS(ci).id;
    do_run = RUN_ALL || eval(['RUN_CONFIG_' id]);

    % Resolve budget for this config
    pop_ci  = CONFIG_POP.(id);  if pop_ci  == 0, pop_ci  = POP_SIZE;   end
    gen_ci  = CONFIG_GEN.(id);  if gen_ci  == 0, gen_ci  = ITERATIONS; end
    runs_ci = CONFIG_RUNS.(id); if runs_ci == 0, runs_ci = N_RUNS;     end
    runs_per_cfg(ci) = runs_ci;

    if do_run
        % ----- FRESH RUN -----
        fprintf('\n=== Config %s: %s | pop=%d gen=%d runs=%d ===\n', ...
            id, CONFIGS(ci).desc, pop_ci, gen_ci, runs_ci);

        use_improved  = CONFIGS(ci).use_improved;
        nsga2_variant = CONFIGS(ci).nsga2_variant;
        init_worker_globals(info.routePath, info.rollingstockPath, cfg.sim, use_improved);
        wait(parfevalOnAll(@init_worker_globals, 0, info.routePath, info.rollingstockPath, ...
            cfg.sim, use_improved));
        dim_all(ci) = dimension;
        pop_size    = pop_ci;
        iterations  = gen_ci;
        wait(parfevalOnAll(@push_worker_globals, 0, pop_ci, gen_ci, SEARCH_TLIM, ...
            CONFIGS(ci).nsga2_variant));

        pops_ci  = cell(runs_ci,1);
        times_ci = zeros(runs_ci,1);
        vp_local = vel_profile;

        parfor (run_id = 1:runs_ci, N_WORKERS)
            rng(run_id, 'twister');
            t0_ = tic;
            [p,~,~] = nsga2_main(vp_local);
            times_ci(run_id) = toc(t0_);
            pops_ci{run_id}  = p;
        end

        for run_id = 1:runs_ci
            pops_all{ci, run_id} = pops_ci{run_id};
            t_mat(ci, run_id)    = times_ci(run_id);
        end
        runs_per_cfg(ci) = runs_ci;
        pop_cfg(ci)      = pop_ci;
        gen_cfg(ci)      = gen_ci;
        fprintf('  Done | dim=%d | mean=%.1fs\n', dim_all(ci), mean(times_ci));

    else
        % ----- LOAD FROM CACHE -----
        fprintf('\n=== Config %s: SKIPPED — loading from %s ===\n', id, mat_cache);
        if exist(mat_cache,'file') ~= 2
            fprintf('  WARNING: cache file not found, skipping Config %s\n', id);
            continue;
        end
        cache = load(mat_cache, 'results');
        rows  = find(strcmp({cache.results.config_id}, id));
        if isempty(rows)
            fprintf('  WARNING: Config %s not found in cache, skipping\n', id);
            continue;
        end
        n_cached = numel(rows);
        runs_per_cfg(ci) = min(runs_ci, n_cached);
        fprintf('  Loaded %d cached runs for Config %s\n', n_cached, id);
        for ri = 1:min(runs_ci, n_cached)
            r = cache.results(rows(ri));
            dim_all(ci) = r.dim;
            t_mat(ci, ri) = r.runtime;
            % Reconstruct pop matrix [AllX | AllF | rank(0) | sde(0)]
            if ~isempty(r.AllF)
                npts = size(r.AllF, 1);
                AX   = r.AllX;
                if isempty(AX) || size(AX,1) ~= npts, AX = zeros(npts, r.dim); end
                pops_all{ci, ri} = [AX, r.AllF, zeros(npts, 2)];
            end
        end
        r0 = cache.results(rows(1));
        pop_cfg(ci) = r0.pop_size;
        gen_cfg(ci) = r0.iterations;
        fprintf('  dim=%d | mean runtime=%.1fs\n', dim_all(ci), mean(t_mat(ci,1:runs_per_cfg(ci)),'omitnan'));
    end
end

%% ===== PHASE 2: UNION PARETO PER CONFIG + GLOBAL REFERENCE =====
% Persis seperti run_benchmark_compare.m:
%   - union semua run per config → kurva Pareto terbaik tiap config
%   - gabung semua → global reference untuk HV/IGD yang comparable
fprintf('\n--- Building union Pareto fronts ---\n');
union_PF  = cell(n_cfg, 1);
allFronts = zeros(0,2);

for ci = 1:n_cfg
    d = dim_all(ci);
    F_union = zeros(0,2);
    for run_id = 1:runs_per_cfg(ci)
        pop_r = pops_all{ci, run_id};
        if isempty(pop_r) || size(pop_r,2) < d+2, continue; end
        F = double(pop_r(:, d+1:d+2));
        valid = F(:,1)<1e5 & F(:,2)<1e5;
        if any(valid), F_union = [F_union; nd_filter(F(valid,:))]; end %#ok<AGROW>
    end
    union_PF{ci} = nd_filter(F_union);
    if ~isempty(union_PF{ci}), allFronts = [allFronts; union_PF{ci}]; end %#ok<AGROW>
    fprintf('  Config %s: union PF = %d pts\n', CONFIGS(ci).id, size(union_PF{ci},1));
end

F_ref_global = nd_filter(allFronts);
fprintf('Global reference front: %d pts\n', size(F_ref_global,1));

mins_r = min(F_ref_global);
maxs_r = max(F_ref_global);
rng_r  = maxs_r - mins_r + eps;
ref_pt = maxs_r + 0.10*(maxs_r - mins_r + eps);

%% ===== PHASE 3: METRICS PER RUN VS GLOBAL REFERENCE =====
hv_mat  = nan(n_cfg, max_runs);
igd_mat = nan(n_cfg, max_runs);
fr_mat  = nan(n_cfg, max_runs);
sv_mat  = nan(n_cfg, max_runs);
E_mat   = nan(n_cfg, max_runs);
all_rows = {};

for ci = 1:n_cfg
    d = dim_all(ci);
    for run_id = 1:runs_per_cfg(ci)
        pop_r = pops_all{ci, run_id};
        if isempty(pop_r) || size(pop_r,2) < d+2
            all_rows{end+1} = {CONFIGS(ci).id, run_id, NaN,NaN,NaN,NaN,NaN,t_mat(ci,run_id)}; %#ok<AGROW>
            continue;
        end
        F_all = double(pop_r(:, d+1:d+2));
        valid = F_all(:,1) < 1e5;
        F_nd  = nd_filter(F_all(valid,:));

        if ~isempty(F_nd) && ~isempty(F_ref_global)
            Fn   = (F_nd - mins_r) ./ rng_r;
            refn = (ref_pt - mins_r) ./ rng_r;
            Rn   = (F_ref_global - mins_r) ./ rng_r;
            hv_mat(ci,run_id)  = hv2d_local(Fn, refn);
            igd_mat(ci,run_id) = igd_local(Rn, Fn);
        end

        fr_mat(ci,run_id) = sum(F_all(valid,1) <= T_TARGET) / max(sum(valid),1);

        E_best = NaN;
        if ~isempty(F_nd)
            feas = F_nd(F_nd(:,1) <= T_TARGET, :);
            if ~isempty(feas), E_best = min(feas(:,2)); end
        end
        E_mat(ci,run_id)  = E_best;
        sv_mat(ci,run_id) = energy_saving_percent(E_best, E_BASELINE);

        all_rows{end+1} = {CONFIGS(ci).id, run_id, hv_mat(ci,run_id), igd_mat(ci,run_id), ...
            fr_mat(ci,run_id), E_best, sv_mat(ci,run_id), t_mat(ci,run_id)}; %#ok<AGROW>
    end

    fprintf('Config %s [%s]\n  HV=%.4f±%.4f | IGD=%.4f | Feas=%.3f | EnSav=%.1f%% | T=%.1fs\n', ...
        CONFIGS(ci).id, CONFIGS(ci).desc, mean(hv_mat(ci,:),'omitnan'), std(hv_mat(ci,:),'omitnan'), ...
        mean(igd_mat(ci,:),'omitnan'), mean(fr_mat(ci,:),'omitnan'), ...
        mean(sv_mat(ci,:),'omitnan'), mean(t_mat(ci,:),'omitnan'));
end

%% ===== SAVE CSV =====
T_csv = cell2table(vertcat(all_rows{:}), ...
    'VariableNames',{'config','run_id','hv','igd','feasible','E_best_kWh','energy_saving_pct','runtime_s'});
writetable(T_csv, fullfile(OUT_DIR,'E1_results.csv'));
fprintf('\nSaved: %s\n', fullfile(OUT_DIR,'E1_results.csv'));

%% ===== SAVE .MAT — kompatibel dengan analyze_cccr_benchmark_full.m =====
% Build full results struct (format sama seperti run_benchmark_compare.m)
tmpl_r = struct('strategy','','solver','','seed',0,'pop_size',0,'dim',0, ...
    'iterations',0,'runtime',NaN,'AllX',zeros(0,0),'AllF',zeros(0,2), ...
    'X',zeros(0,0),'F',zeros(0,2),'HV',NaN,'IGD',NaN,'Nf1',0, ...
    'nsga2_variant','','config_id','','config_desc','', ...
    'use_improved',false,'route_file','','rs_file','');

all_results = repmat(tmpl_r, 0, 1);
for ci = 1:n_cfg
    d = dim_all(ci);
    for run_id = 1:runs_per_cfg(ci)
        r = tmpl_r;
        r.solver        = 'nsga2';
        r.seed          = run_id;
        r.pop_size      = pop_cfg(ci);
        r.dim           = d;
        r.iterations    = gen_cfg(ci);
        r.runtime       = t_mat(ci, run_id);
        r.HV            = hv_mat(ci, run_id);
        r.IGD           = igd_mat(ci, run_id);
        r.nsga2_variant = CONFIGS(ci).nsga2_variant;
        r.config_id     = CONFIGS(ci).id;
        r.config_desc   = CONFIGS(ci).desc;
        r.use_improved  = CONFIGS(ci).use_improved;
        r.strategy      = ternary(CONFIGS(ci).use_improved, 'improved', 'base');
        r.route_file    = ROUTE_IS04;
        r.rs_file       = RS_FILE;

        pop_r = pops_all{ci, run_id};
        if ~isempty(pop_r) && size(pop_r,2) >= d+2
            r.AllX    = double(pop_r(:, 1:d));
            r.AllF    = double(pop_r(:, d+1:d+2));
            valid     = r.AllF(:,1)<1e5 & r.AllF(:,2)<1e5;
            X_valid   = r.AllX(valid,:);
            F_valid   = r.AllF(valid,:);
            [r.F, nd_idx] = nd_filter_idx(F_valid);  % Pareto [T E] + row indices
            r.X       = X_valid(nd_idx,:);            % X for each Pareto point
            r.Nf1     = size(r.F, 1);
        end
        all_results(end+1,1) = r; %#ok<AGROW>
    end
end

% Simpan semua 4 config (A, B, C, D) dalam satu file
results = all_results;
save(fullfile(OUT_DIR,'E1_all_configs.mat'), 'results');
fprintf('Saved: %s\n', fullfile(OUT_DIR,'E1_all_configs.mat'));

%% ===== SUMMARY + WILCOXON =====
fid = fopen(fullfile(OUT_DIR,'E1_summary.txt'),'w');
fprintf(fid,'Experiment E1 — Ablation Study\n');
fprintf(fid,'Segment: %s | T_target=%d s | E_baseline=%.4f kWh\n', SEG_LABEL, T_TARGET, E_BASELINE);
fprintf(fid,'Pop=%d | Gen=%d | Runs=%d | Workers=%d\n\n', POP_SIZE, ITERATIONS, N_RUNS, N_WORKERS);
fprintf(fid,'%-45s  %8s  %8s  %8s  %10s  %10s\n', ...
    'Config','HV_mean','HV_std','IGD_mean','Feas_mean','EnSav_pct');

for ci = 1:n_cfg
    line = sprintf('%-45s  %8.4f  %8.4f  %8.4f  %10.3f  %10.2f', ...
        [CONFIGS(ci).id ': ' CONFIGS(ci).desc], ...
        mean(hv_mat(ci,:),'omitnan'), std(hv_mat(ci,:),'omitnan'), ...
        mean(igd_mat(ci,:),'omitnan'), mean(fr_mat(ci,:),'omitnan'), ...
        mean(sv_mat(ci,:),'omitnan'));
    disp(line); fprintf(fid,'%s\n',line);
end

HV_A = mean(hv_mat(1,:),'omitnan'); HV_B = mean(hv_mat(2,:),'omitnan');
HV_C = mean(hv_mat(3,:),'omitnan'); HV_D = mean(hv_mat(4,:),'omitnan');
fprintf(fid,'\nRelative HV contributions:\n');
fprintf(fid,'  Strategy : %+.2f%%\n', (HV_C-HV_A)/(HV_A+eps)*100);
fprintf(fid,'  Solver   : %+.2f%%\n', (HV_B-HV_A)/(HV_A+eps)*100);
fprintf(fid,'  Interact : %+.2f%%\n', ((HV_D-HV_A)-(HV_C-HV_A)-(HV_B-HV_A))/(HV_A+eps)*100);

fprintf(fid,'\nWilcoxon signed-rank test (HV, paired — same seeds):\n');
pairs = {1,4; 1,2; 1,3};
for pi = 1:size(pairs,1)
    c1=pairs{pi,1}; c2=pairs{pi,2};
    x=hv_mat(c1,:)'; y=hv_mat(c2,:)';
    ok=~isnan(x)&~isnan(y);
    if sum(ok)>=5 && exist('signrank','file')==2
        p_val=signrank(x(ok),y(ok));
    else
        p_val=wilcoxon_srtest(x(ok),y(ok));
    end
    line=sprintf('  [%s] %s\n    vs [%s] %s\n    p = %.4f%s', ...
        CONFIGS(c1).id, CONFIGS(c1).desc, CONFIGS(c2).id, CONFIGS(c2).desc, ...
        p_val, ternary(p_val<0.05,'  *',''));
    disp(line); fprintf(fid,'%s\n',line);
end

% Energy saving per comparison pair
fprintf(fid,'\nEnergy saving comparison:\n');
pairs_lbl = {'A vs C (Strategy effect)','A vs B (Solver effect)','A vs D (Combined effect)'};
pairs_idx  = {[1 3],[1 2],[1 4]};
for pi = 1:3
    ca=pairs_idx{pi}(1); cb=pairs_idx{pi}(2);
    sav_a=mean(sv_mat(ca,:),'omitnan');
    sav_b=mean(sv_mat(cb,:),'omitnan');
    fprintf(fid,'  %s: %.2f%% vs %.2f%% (diff=%.2f%%)\n', ...
        pairs_lbl{pi}, sav_a, sav_b, sav_b-sav_a);
end
fclose(fid);

%% ===== PLOT 1: PARETO FRONT OVERLAY =====
% Semua 4 config dalam 1 grafik — seperti plot_pareto_overlay di run_benchmark_compare.m
colors_cfg = [0.55 0.55 0.55; 0.2 0.5 0.9; 0.9 0.55 0.1; 0.15 0.72 0.15];
lw_cfg     = [1.2, 1.5, 1.5, 2.5];

figure(1001); clf; hold on;
for ci = 1:n_cfg
    F = union_PF{ci};
    if isempty(F), continue; end
    [~,ord] = sort(F(:,1));
    plot(F(ord,1), F(ord,2), '-o', 'Color', colors_cfg(ci,:), ...
        'LineWidth', lw_cfg(ci), 'MarkerSize', 5, ...
        'DisplayName', [CONFIGS(ci).id ': ' CONFIGS(ci).desc]);
end
xline(T_TARGET,'--k','T_{target}','LabelVerticalAlignment','bottom','LineWidth',1.5);
xlim([0, PLOT_TLIM]);
xlabel('Running time (s)'); ylabel('Energy (kWh)');
title(sprintf('E1: Pareto Front Overlay — %s | T_{target}=%ds', SEG_LABEL, T_TARGET));
legend('Location','northeast','FontSize',9); grid on;
saveas(gcf, fullfile(OUT_DIR,'E1_pareto_overlay.png'));

%% ===== PLOT 2: ENERGY ENVELOPE COMPARISON =====
% Seperti plot_energy_compare di run_benchmark_compare.m
% Panel kiri: effect of strategy (A vs C)
% Panel kanan: effect of solver (A vs B)
figure(1002); clf;
subplot(1,2,1);
plot_energy_envelope(union_PF{1}, union_PF{3}, T_TARGET, PLOT_TLIM, ...
    ['A: ' CONFIGS(1).desc], ['C: ' CONFIGS(3).desc], ...
    'Effect of Driving Strategy (A vs C)');

subplot(1,2,2);
plot_energy_envelope(union_PF{1}, union_PF{2}, T_TARGET, PLOT_TLIM, ...
    ['A: ' CONFIGS(1).desc], ['B: ' CONFIGS(2).desc], ...
    'Effect of Solver (A vs B)');
sgtitle(sprintf('E1: Energy Saving Comparison — %s', SEG_LABEL));
saveas(gcf, fullfile(OUT_DIR,'E1_energy_comparison.png'));

% Grafik terpisah: A vs D (combined/proposed vs baseline)
figure(1003); clf;
plot_energy_envelope(union_PF{1}, union_PF{4}, T_TARGET, PLOT_TLIM, ...
    ['A: ' CONFIGS(1).desc], ['D: ' CONFIGS(4).desc], ...
    sprintf('E1: Combined Effect (A vs D) — %s', SEG_LABEL));
saveas(gcf, fullfile(OUT_DIR,'E1_energy_combined.png'));

%% ===== PLOT 3: HV BAR CHART =====
figure(1004); clf;
hv_means = mean(hv_mat, 2, 'omitnan');
hv_stds  = std(hv_mat, 0, 2, 'omitnan');
b = bar(hv_means, 0.6, 'FaceColor','flat');
colors_bar = [0.7 0.7 0.7; 0.5 0.7 0.9; 0.9 0.7 0.5; 0.4 0.8 0.4];
for ci=1:n_cfg, b.CData(ci,:) = colors_bar(ci,:); end
hold on;
errorbar(1:n_cfg, hv_means, hv_stds, 'k.', 'LineWidth',1.5, 'CapSize',8);
set(gca,'XTickLabel',{CONFIGS.id},'XTick',1:n_cfg,'FontSize',11);
xlabel('Configuration'); ylabel('Hypervolume (HV)');
title(sprintf('E1: Mean HV ± std | %s | T_{target}=%ds | %d runs', SEG_LABEL, T_TARGET, N_RUNS));
legend({'Mean HV','±1 std'},'Location','northwest');
grid on; box on;
saveas(gcf, fullfile(OUT_DIR,'E1_hv_bar.png'));

fprintf('\nE1 complete → %s\n', OUT_DIR);

%% ===== HELPERS =====
function plot_energy_envelope(F_base, F_comp, T_tgt, T_lim, lbl_base, lbl_comp, ttl)
    env_b = make_envelope(F_base);
    env_c = make_envelope(F_comp);
    if ~isempty(env_b.T)
        plot(env_b.T, env_b.E, '-o', 'LineWidth',1.5, 'DisplayName', lbl_base);
    end
    if ~isempty(env_c.T)
        plot(env_c.T, env_c.E, '-s', 'LineWidth',2.0, 'DisplayName', lbl_comp);
    end
    if ~isempty(env_b.T) && ~isempty(env_c.T)
        Tmin=max(min(env_b.T),min(env_c.T)); Tmax=min(max(env_b.T),max(env_c.T));
        if Tmin < Tmax
            Tg=linspace(Tmin,Tmax,60)';
            Eb=interp1(env_b.T,env_b.E,Tg,'pchip','extrap');
            Ec=interp1(env_c.T,env_c.E,Tg,'pchip','extrap');
            sav=100*(Eb-Ec)./max(Eb,eps);
            fprintf('%s: mean saving = %.2f%%\n', ttl, mean(sav,'omitnan'));
        end
    end
    xline(T_tgt,'--k','T_{target}','LabelVerticalAlignment','bottom');
    xlim([0, T_lim]);
    xlabel('Running time (s)'); ylabel('Best energy (kWh)');
    title(ttl); legend('Location','best','FontSize',9); grid on;
end

function env = make_envelope(F)
    if isempty(F), env=struct('T',zeros(0,1),'E',zeros(0,1)); return; end
    [Ts,ord]=sort(F(:,1)); Es=F(ord,2);
    Eb=zeros(size(Es)); m=inf;
    for i=1:numel(Es)
        if Es(i)<m, m=Es(i); end; Eb(i)=m;
    end
    [Tu,ia]=unique(Ts,'stable');
    env=struct('T',Tu,'E',Eb(ia));
end

function F_nd = nd_filter(F)
    if isempty(F), F_nd=zeros(0,2); return; end
    n=size(F,1); k=true(n,1);
    for i=1:n
        if ~k(i),continue;end
        for j=1:n
            if i==j||~k(j),continue;end
            if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)),k(i)=false;break;end
        end
    end
    F_nd=F(k,:); [~,o]=sort(F_nd(:,1)); F_nd=F_nd(o,:);
end

function hv = hv2d_local(F, ref)
    if isempty(F), hv=0; return; end
    F=nd_filter(F); F(:,2)=cummin(F(:,2));
    hv=0; pE=ref(2);
    for i=1:size(F,1)
        if F(i,2)>pE, continue; end
        hv=hv+(ref(1)-F(i,1))*(pE-F(i,2)); pE=F(i,2);
    end
end

function val = igd_local(R, A)
    if isempty(R)||isempty(A), val=NaN; return; end
    d=zeros(size(R,1),1);
    for i=1:size(R,1)
        dr=A-R(i,:); d(i)=min(sqrt(sum(dr.^2,2)));
    end
    val=mean(d);
end

function v = ternary(c,a,b)
    if c, v=a; else, v=b; end
end

function [F_nd, idx] = nd_filter_idx(F)
% Like nd_filter but also returns row indices into F so caller can extract X.
    if isempty(F), F_nd=zeros(0,2); idx=zeros(0,1); return; end
    n=size(F,1); k=true(n,1);
    for i=1:n
        if ~k(i),continue;end
        for j=1:n
            if i==j||~k(j),continue;end
            if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)),k(i)=false;break;end
        end
    end
    idx_raw = find(k);
    [~,o]   = sort(F(idx_raw,1));
    idx     = idx_raw(o);
    F_nd    = F(idx,:);
end

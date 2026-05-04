%% run_E2_tight_sweep.m — Experiment E2: Tight Journey Time Sweep (G1)
%
% Config D (improved CC_CR + BRL-SDE-NSGA-II) across 8 slack levels.
% Segments: IS04 (1642 m) and IS08 (3778 m).
% Slack   : [1, 2, 5, 8.4, 10, 15, 20, 30] %
% T_target: T_min × (1 + slack/100)
% Runs    : 30 per (segment, slack) combination.
% Outputs : E2_results.csv, E2_summary.txt, two figures

clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% ===== SETTINGS =====
RS_FILE    = 'rollingstock_Guangzhou_L7.m';
SLACK_VALS = [1, 2, 5, 8.4, 10, 15, 20, 30];
POP_SIZE   = 200;
ITERATIONS = 300;
FEAS_CRIT  = 0.95;

% Jika dipanggil dari main_experiments.m, gunakan ACTIVE_SEG dari workspace
if exist('ACTIVE_SEG','var') && ~isempty(ACTIVE_SEG)
    SEGMENTS = struct('name', {ACTIVE_SEG.name}, 'file', {ACTIVE_SEG.file});
    N_RUNS   = ACTIVE_NRUNS;
    OUT_DIR  = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', ACTIVE_SEG.name);
else
    SEGMENTS = struct( ...
        'name', {'IS04','IS08'}, ...
        'file', {'Guangzhou_Line7_IS04_5.200-6.842km.mat', ...
                 'Guangzhou_Line7_IS08_13.729-17.507km.mat'});
    N_RUNS   = 30;
    OUT_DIR  = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results');
end
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

%% ===== GLOBALS =====
global use_improved nsga2_variant pop_size iterations dimension
global vel_profile time_obj_max parallel_use show_progress driving_strategy

use_improved     = true;
nsga2_variant    = 'rl_sde';
show_progress    = false;
parallel_use     = true;
driving_strategy = "CC_CR";

%% ===== LOOP OVER SEGMENTS =====
all_rows = {};
n_seg    = numel(SEGMENTS);
n_sl     = numel(SLACK_VALS);

mean_hv   = nan(n_seg, n_sl);
mean_feas = nan(n_seg, n_sl);
mean_esav = nan(n_seg, n_sl);
mean_igd  = nan(n_seg, n_sl);
mean_fsz  = nan(n_seg, n_sl);

for si = 1:n_seg
    seg_name = SEGMENTS(si).name;
    seg_file = SEGMENTS(si).file;
    fprintf('\n========== Segment: %s ==========\n', seg_name);

    cfg = struct('route_file', seg_file, 'rollingstock_file', RS_FILE, ...
        'driving_strategy', "CC_CR", 'pop_size', POP_SIZE, 'iterations', ITERATIONS, ...
        'use_improved', true, 'nsga2_variant', 'rl_sde', 'parallel_use', true, ...
        'sim', struct());
    info = setup_project(cfg);

    setup_parallel_pool(cfg, info);
    pool = gcp('nocreate');
    if ~isempty(pool)
        wait(parfevalOnAll(@init_worker_globals, 0, info.routePath, info.rollingstockPath, ...
            cfg.sim, true));
    end

    [T_min, E_base] = flat_out_baseline_rk4(info.routePath, info.rollingstockPath);
    fprintf('  T_min=%.2f s | E_baseline=%.4f kWh\n', T_min, E_base);

    for sli = 1:n_sl
        slack    = SLACK_VALS(sli);
        T_target = T_min * (1 + slack/100);
        time_obj_max = T_target * 1.50;   % lebih longgar — mode santai tetap bisa ditemukan
        pop_size     = POP_SIZE;
        iterations   = ITERATIONS;

        fprintf('\n  -- Slack=%.1f%% | T_target=%.2f s --\n', slack, T_target);

        % Push all three globals to workers
        if ~isempty(pool)
            wait(parfevalOnAll(@push_worker_globals, 0, POP_SIZE, ITERATIONS, T_target*1.50));
        end

        hv_v  = nan(N_RUNS,1); igd_v = nan(N_RUNS,1);
        fr_v  = nan(N_RUNS,1); es_v  = nan(N_RUNS,1);
        rt_v  = nan(N_RUNS,1); fs_v  = nan(N_RUNS,1);
        pops  = cell(N_RUNS,1); tims  = zeros(N_RUNS,1);
        all_F_ref = zeros(0,2);

        parfor (run_id = 1:N_RUNS, 4)
            rng(run_id, 'twister');
            t0_ = tic;
            [p,~,~] = nsga2_main(vel_profile);
            tims(run_id) = toc(t0_);
            pops{run_id} = p;
        end

        for run_id = 1:N_RUNS
            F_a = double(pops{run_id}(:, dimension+1:dimension+2));
            valid = F_a(:,1) < 1e5;
            if any(valid), all_F_ref=[all_F_ref; nd_filter(F_a(valid,:))]; end %#ok<AGROW>
        end
        F_ref = nd_filter(all_F_ref);

        for run_id = 1:N_RUNS
            F_all = double(pops{run_id}(:, dimension+1:dimension+2));
            valid = F_all(:,1) < 1e5;
            F_nd  = nd_filter(F_all(valid,:));

            [hv_v(run_id), igd_v(run_id), ~, fr_v(run_id)] = ...
                compute_metrics(F_all, F_ref, T_target);
            rt_v(run_id) = tims(run_id);
            fs_v(run_id) = size(F_nd,1);

            E_best = NaN;
            if ~isempty(F_nd)
                feas_nd = F_nd(F_nd(:,1) <= T_target, :);
                if ~isempty(feas_nd), E_best = min(feas_nd(:,2)); end
            end
            es_v(run_id) = energy_saving_percent(E_best, E_base);

            all_rows{end+1} = {seg_name, slack, run_id, ...
                hv_v(run_id), igd_v(run_id), fr_v(run_id), ...
                es_v(run_id), fs_v(run_id), rt_v(run_id)}; %#ok<AGROW>
        end

        mean_hv(si,sli)   = mean(hv_v,  'omitnan');
        mean_feas(si,sli) = mean(fr_v,  'omitnan');
        mean_esav(si,sli) = mean(es_v,  'omitnan');
        mean_igd(si,sli)  = mean(igd_v, 'omitnan');
        mean_fsz(si,sli)  = mean(fs_v,  'omitnan');

        fprintf('    HV=%.4f | Feas=%.3f | EnSav=%.1f%% | FrontSz=%.1f\n', ...
            mean_hv(si,sli), mean_feas(si,sli), mean_esav(si,sli), mean_fsz(si,sli));
    end
end

%% ===== SAVE CSV =====
T_csv = cell2table(vertcat(all_rows{:}), ...
    'VariableNames',{'segment','slack_pct','run_id','hv','igd','feasible', ...
                     'energy_saving_pct','front_size','runtime_s'});
writetable(T_csv, fullfile(OUT_DIR,'E2_results.csv'));
fprintf('\nSaved: %s\n', fullfile(OUT_DIR,'E2_results.csv'));

%% ===== SUMMARY =====
fid = fopen(fullfile(OUT_DIR,'E2_summary.txt'),'w');
fprintf(fid,'Experiment E2 — Tight Journey Time Sweep\n');
fprintf(fid,'Config D | Slack=[%.4g] %%\n\n', SLACK_VALS);
fprintf(fid,'%-10s', 'Segment');
for sli=1:n_sl, fprintf(fid,'%9.4g%%', SLACK_VALS(sli)); end
fprintf(fid,'\n');
for si=1:n_seg
    fprintf(fid,'%-10s', SEGMENTS(si).name);
    for sli=1:n_sl, fprintf(fid,'%9.2f', mean_esav(si,sli)); end
    fprintf(fid,' (EnSav%%)\n');
end
fprintf(fid,'\nCritical slack (feasibility drops below %.0f%%):\n', FEAS_CRIT*100);
for si=1:n_seg
    crit = find(mean_feas(si,:) < FEAS_CRIT, 1, 'first');
    if ~isempty(crit)
        fprintf(fid,'  %s: critical at slack=%.4g%%  (feas=%.3f)\n', ...
            SEGMENTS(si).name, SLACK_VALS(crit), mean_feas(si,crit));
    else
        fprintf(fid,'  %s: feasibility >= %.0f%% across all tested slack levels\n', ...
            SEGMENTS(si).name, FEAS_CRIT*100);
    end
end
fclose(fid);

fprintf('\n--- Energy Saving by Slack ---\n');
fprintf('%-6s', 'Seg');
for sli=1:n_sl, fprintf('  %5.4g%%', SLACK_VALS(sli)); end; fprintf('\n');
for si=1:n_seg
    fprintf('%-6s', SEGMENTS(si).name);
    for sli=1:n_sl, fprintf('  %5.1f', mean_esav(si,sli)); end; fprintf('\n');
end

%% ===== PLOTS =====
colors = lines(n_seg);

figure(2001); clf; hold on;
for si=1:n_seg
    plot(SLACK_VALS, mean_esav(si,:), '-o', 'Color',colors(si,:), ...
        'LineWidth',1.5, 'DisplayName', SEGMENTS(si).name);
end
xlabel('Slack (%)'); ylabel('Mean Energy Saving (%)');
title('E2: Energy Saving vs Slack — Config D (improved + BRL-SDE)');
legend('Location','northwest'); grid on; xlim([0 32]);
saveas(gcf, fullfile(OUT_DIR,'E2_energy_saving_vs_slack.png'));

figure(2002); clf; hold on;
for si=1:n_seg
    plot(SLACK_VALS, mean_feas(si,:), '-s', 'Color',colors(si,:), ...
        'LineWidth',1.5, 'DisplayName', SEGMENTS(si).name);
end
yline(FEAS_CRIT,'--r', sprintf('Critical = %.0f%%', FEAS_CRIT*100), 'LineWidth',1.5);
xlabel('Slack (%)'); ylabel('Mean Feasibility Rate');
title('E2: Feasibility Rate vs Slack — Config D');
legend('Location','southeast'); grid on;
ylim([0 1.05]); xlim([0 32]);
for si=1:n_seg
    crit=find(mean_feas(si,:)<FEAS_CRIT,1,'first');
    if ~isempty(crit)
        plot(SLACK_VALS(crit), mean_feas(si,crit), 'rv', 'MarkerSize',10, ...
            'HandleVisibility','off');
    end
end
saveas(gcf, fullfile(OUT_DIR,'E2_feasibility_vs_slack.png'));

fprintf('\nE2 complete → %s\n', OUT_DIR);

%% ===== HELPERS =====
function F_nd = nd_filter(F)
    if isempty(F), F_nd=zeros(0,2); return; end
    n=size(F,1); k=true(n,1);
    for i=1:n
        if ~k(i), continue; end
        for j=1:n
            if i==j||~k(j), continue; end
            if all(F(j,:)<=F(i,:))&&any(F(j,:)<F(i,:)), k(i)=false; break; end
        end
    end
    F_nd=F(k,:); [~,o]=sort(F_nd(:,1)); F_nd=F_nd(o,:);
end

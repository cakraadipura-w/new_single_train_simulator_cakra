function analyze_E1_thresholds(seg_name, T_thresholds)
%% analyze_E1_thresholds — Re-analyze E1 Pareto fronts at various T_threshold
%
% Tujuan: Tanpa run ulang NSGA-II, cek apakah Config D menang vs A/B/C
%         pada berbagai nilai T_threshold (dari longgar ke ketat).
%
% Usage:
%   analyze_E1_thresholds                     % default: IS02, auto threshold
%   analyze_E1_thresholds('IS04')             % segment lain
%   analyze_E1_thresholds('IS02', [110 120 130 140 170])   % threshold custom

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

if nargin < 1 || isempty(seg_name), seg_name = 'IS02'; end

%% --- Load MAT file ---
mat_file = fullfile(pwd, 'experiment_results', seg_name, 'E1_all_configs.mat');
if ~exist(mat_file, 'file')
    mat_file = fullfile(pwd, 'experiment_results', seg_name, 'E1_A_vs_D.mat');
end
assert(exist(mat_file,'file')==2, 'File tidak ditemukan: %s', mat_file);

data    = load(mat_file);
results = data.results;
fprintf('Loaded: %s  (%d result entries)\n', mat_file, numel(results));

%% --- Estimasi T_min dari data Pareto ---
all_T = [];
for i = 1:numel(results)
    if ~isempty(results(i).F)
        all_T = [all_T; results(i).F(:,1)]; %#ok<AGROW>
    end
end
T_min_est = min(all_T);
fprintf('T_min (estimated from Pareto data) = %.2f s\n\n', T_min_est);

%% --- Default thresholds: dari ketat ke longgar ---
if nargin < 2 || isempty(T_thresholds)
    slack_pcts   = [5, 10, 15, 20, 25, 30, 40, 50, 70];
    T_thresholds = T_min_est * (1 + slack_pcts/100);
end

%% --- Ambil config unik ---
cfg_ids  = unique({results.config_id}, 'stable');
n_cfg    = numel(cfg_ids);
n_thr    = numel(T_thresholds);
colors   = lines(n_cfg);

%% --- Hitung E_best, std, feas_rate per config per threshold ---
E_best_mean = nan(n_cfg, n_thr);
E_best_std  = nan(n_cfg, n_thr);
feas_mean   = nan(n_cfg, n_thr);
n_runs_cfg  = zeros(n_cfg, 1);

for ci = 1:n_cfg
    runs = find(strcmp({results.config_id}, cfg_ids{ci}));
    n_runs_cfg(ci) = numel(runs);

    for ti = 1:n_thr
        T_thr = T_thresholds(ti);
        E_run = nan(numel(runs), 1);
        f_run = nan(numel(runs), 1);

        for ri = 1:numel(runs)
            F = results(runs(ri)).F;   % [T, E] Pareto front
            if isempty(F), continue; end
            feas_pts  = F(F(:,1) <= T_thr, :);
            total_pts = size(F, 1);
            f_run(ri) = size(feas_pts,1) / total_pts;
            if ~isempty(feas_pts)
                E_run(ri) = min(feas_pts(:,2));
            end
        end
        E_best_mean(ci,ti) = mean(E_run,  'omitnan');
        E_best_std(ci,ti)  = std(E_run, 0,'omitnan');
        feas_mean(ci,ti)   = mean(f_run,  'omitnan');
    end
end

%% --- Print tabel per threshold ---
fprintf('%-10s  %-8s  %-10s  %-8s  %-8s\n', ...
    'T_thr(s)', 'Slack%', 'Config', 'E_best', 'Feas%');
fprintf('%s\n', repmat('-',1,55));

for ti = 1:n_thr
    T_thr     = T_thresholds(ti);
    slack_pct = (T_thr - T_min_est) / T_min_est * 100;
    fprintf('\nT_threshold = %.1f s  (slack %.1f%% diatas T_min)\n', T_thr, slack_pct);
    fprintf('  %-8s  %10s ± %-8s  %8s\n','Config','E_best','std','Feas%');

    [~, ord] = sort(E_best_mean(:,ti));
    for ri = 1:n_cfg
        ci = ord(ri);
        mark = '';
        if ri == 1, mark = '  ← BEST E'; end
        fprintf('  %-8s  %10.4f ± %-8.4f  %7.1f%%%s\n', ...
            cfg_ids{ci}, E_best_mean(ci,ti), E_best_std(ci,ti), ...
            feas_mean(ci,ti)*100, mark);
    end
end

%% --- Plot 1: E_best vs T_threshold ---
figure('Name', sprintf('E1 Re-analysis: E_best vs Threshold — %s', seg_name), ...
    'NumberTitle','off', 'Position',[100 100 800 450]);
clf; hold on;

for ci = 1:n_cfg
    lw = 2.5 + (strcmp(cfg_ids{ci},'D')*0.5);   % D lebih tebal
    plot(T_thresholds, E_best_mean(ci,:), '-o', ...
        'Color', colors(ci,:), 'LineWidth', lw, ...
        'MarkerFaceColor', colors(ci,:), ...
        'DisplayName', sprintf('Config %s', cfg_ids{ci}));
    errorbar(T_thresholds, E_best_mean(ci,:), E_best_std(ci,:), ...
        'Color', colors(ci,:), 'LineStyle','none', 'HandleVisibility','off');
end

xline(T_min_est*1.20, '--k', 'T_{min}×1.20 (slack 20%)', ...
    'LabelVerticalAlignment','bottom', 'LineWidth',1.5);
xlabel('T_{threshold} (s)  — ambil Pareto points: T ≤ threshold');
ylabel('E_{best} mean (kWh)');
title(sprintf('%s — Config mana yang menang di setiap T threshold?', seg_name));
legend('Location','best'); grid on; box on;

%% --- Plot 2: Feasibility rate vs T_threshold ---
figure('Name', sprintf('E1 Re-analysis: Feasibility — %s', seg_name), ...
    'NumberTitle','off', 'Position',[950 100 800 450]);
clf; hold on;

for ci = 1:n_cfg
    lw = 2.5 + (strcmp(cfg_ids{ci},'D')*0.5);
    plot(T_thresholds, feas_mean(ci,:)*100, '-s', ...
        'Color', colors(ci,:), 'LineWidth', lw, ...
        'MarkerFaceColor', colors(ci,:), ...
        'DisplayName', sprintf('Config %s', cfg_ids{ci}));
end

xline(T_min_est*1.20, '--k', 'T_{min}×1.20', 'LabelVerticalAlignment','bottom');
xlabel('T_{threshold} (s)');
ylabel('% Pareto points yang T ≤ threshold (feas rate)');
title(sprintf('%s — Feasibility rate tiap Config vs T threshold', seg_name));
legend('Location','best'); grid on; box on;

%% --- Plot 3: Pareto fronts overlay semua config ---
figure('Name', sprintf('E1 Pareto Front Overlay — %s', seg_name), ...
    'NumberTitle','off', 'Position',[100 600 900 450]);
clf; hold on;

% aggregate Pareto per config (union of all runs, then nd-filter)
for ci = 1:n_cfg
    runs = find(strcmp({results.config_id}, cfg_ids{ci}));
    F_all_cfg = [];
    for ri = 1:numel(runs)
        F = results(runs(ri)).F;
        if ~isempty(F), F_all_cfg = [F_all_cfg; F]; end %#ok<AGROW>
    end
    F_nd = nd_filter_local(F_all_cfg);
    if ~isempty(F_nd)
        [~, ord] = sort(F_nd(:,1));
        F_nd = F_nd(ord,:);
        lw = 2.5 + (strcmp(cfg_ids{ci},'D')*0.5);
        plot(F_nd(:,1), F_nd(:,2), '-', ...
            'Color', colors(ci,:), 'LineWidth', lw, ...
            'DisplayName', sprintf('Config %s', cfg_ids{ci}));
    end
end

% Tandai T_min*1.20
T_line = T_min_est * 1.20;
xline(T_line, '--k', sprintf('T_{min}×1.20 = %.1f s', T_line), ...
    'LabelVerticalAlignment','bottom', 'LineWidth',1.5);
xlabel('Running time T (s)');
ylabel('Energy E (kWh)');
title(sprintf('%s — Gabungan Pareto front semua runs per Config', seg_name));
legend('Location','best'); grid on; box on;
xlim([T_min_est*0.98, max(T_thresholds)*1.05]);

fprintf('\nDone. 3 figures generated.\n');
end

%% ===== LOCAL HELPER =====
function F_nd = nd_filter_local(F)
    if isempty(F), F_nd = zeros(0,2); return; end
    n = size(F,1); k = true(n,1);
    for i = 1:n
        if ~k(i), continue; end
        for j = 1:n
            if i==j || ~k(j), continue; end
            if all(F(j,:)<=F(i,:)) && any(F(j,:)<F(i,:))
                k(i) = false; break;
            end
        end
    end
    F_nd = F(k,:);
end

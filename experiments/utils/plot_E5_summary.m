function plot_E5_summary(seg_name, T_sched)
%% PLOT_E5_SUMMARY  Re-display E5 results from the saved MAT file only.
%
% Usage:
%   plot_E5_summary                  % default: IS02
%   plot_E5_summary('IS08')
%   plot_E5_summary('IS02', 170)

root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root_dir));

if nargin < 1 || isempty(seg_name), seg_name = 'IS03'; end
if nargin < 2, T_sched = []; end

T_SCHED_LUT = struct('IS01',130,'IS02',170,'IS03',185,'IS04',180, ...
                     'IS05',185,'IS06',220,'IS07',210,'IS08',330);

res_dir  = fullfile(pwd, 'experiment_results', seg_name);
mat_file = fullfile(res_dir, 'E5_all_configs.mat');
assert(exist(mat_file,'file')==2, 'File not found: %s', mat_file);

data = load(mat_file);
assert(isfield(data, 'results'), 'Variable "results" not found in %s', mat_file);

results = data.results;
meta = local_get_field(data, 'meta', struct());
advisor_request = local_get_field(data, 'advisor_request', struct());
metrics_table = local_get_field(data, 'metrics_table', table());

cfg_ids = {'A','B','C','D','C_LLM','D_LLM'};
n_cfg = numel(cfg_ids);

CLR = struct('A',[0.55 0.55 0.55], 'B',[0.20 0.45 0.85], ...
             'C',[0.95 0.50 0.05], 'D',[0.08 0.78 0.18], ...
             'C_LLM',[0.95 0.70 0.25], 'D_LLM',[0.15 0.55 0.15]);
MRK = struct('A','o', 'B','s', 'C','^', 'D','d', 'C_LLM','v', 'D_LLM','p');
LBL = struct( ...
    'A', 'A: Orig. CC_CR + Orig. NSGA-II', ...
    'B', 'B: Orig. CC_CR + BRL-SDE', ...
    'C', 'C: Improved CC_CR + Orig. NSGA-II', ...
    'D', 'D: Improved CC_CR + BRL-SDE', ...
    'C_LLM', 'C_LLM: C + advisor prior', ...
    'D_LLM', 'D_LLM: D + advisor prior');

segment_name = local_pick_text(seg_name, local_get_field(meta, 'segment_name', ''), local_get_field(advisor_request, 'segment_name', ''));
if isempty(T_sched)
    T_sched = local_get_field(meta, 'T_target_s', NaN);
end
if (isempty(T_sched) || isnan(T_sched)) && isfield(advisor_request, 't_target_s')
    T_sched = advisor_request.t_target_s;
end
if (isempty(T_sched) || isnan(T_sched)) && isfield(T_SCHED_LUT, segment_name)
    T_sched = T_SCHED_LUT.(segment_name);
end

E_baseline = local_get_field(meta, 'E_baseline_kWh', NaN);
if isnan(E_baseline)
    baseline_file = fullfile(res_dir, 'baseline.mat');
    if exist(baseline_file, 'file') == 2
        base_data = load(baseline_file);
        if isfield(base_data, 'E_base_seg')
            E_baseline = base_data.E_base_seg;
        end
    end
end

all_T = [];
for i = 1:numel(results)
    F_here = local_valid_front(results(i));
    if ~isempty(F_here)
        all_T = [all_T; F_here(:,1)]; %#ok<AGROW>
    end
end
assert(~isempty(all_T), 'No valid Pareto data found in %s', mat_file);
T_min_est = min(all_T);

hv_mean   = nan(n_cfg,1);
hv_std    = nan(n_cfg,1);
igd_mean  = nan(n_cfg,1);
igd_std   = nan(n_cfg,1);
feas_mean = nan(n_cfg,1);
E_mean    = nan(n_cfg,1);
runtime_m = nan(n_cfg,1);
ensav_m   = nan(n_cfg,1);
run_count = zeros(n_cfg,1);

for ci = 1:n_cfg
    cfg_runs = results(strcmp({results.config_id}, cfg_ids{ci}));
    run_count(ci) = numel(cfg_runs);
    hv_v = [cfg_runs.HV]';
    igd_v = [cfg_runs.IGD]';
    rt_v = [cfg_runs.runtime]';
    feas_v = nan(numel(cfg_runs),1);
    e_v = nan(numel(cfg_runs),1);
    for ri = 1:numel(cfg_runs)
        F_here = local_valid_front(cfg_runs(ri));
        if isempty(F_here)
            continue;
        end
        feas = F_here(F_here(:,1) <= T_sched, :);
        feas_v(ri) = ~isempty(feas);
        if ~isempty(feas)
            e_v(ri) = min(feas(:,2));
        end
    end

    if ~isempty(metrics_table)
        mask = strcmp(metrics_table.config, cfg_ids{ci});
        if any(mask)
            hv_v = metrics_table.hv(mask);
            igd_v = metrics_table.igd(mask);
            rt_v = metrics_table.runtime_s(mask);
            if ismember('feasible', metrics_table.Properties.VariableNames)
                feas_v = metrics_table.feasible(mask);
            end
            if ismember('E_best_kWh', metrics_table.Properties.VariableNames)
                e_v = metrics_table.E_best_kWh(mask);
            end
        end
    end

    hv_mean(ci)   = mean(hv_v, 'omitnan');
    hv_std(ci)    = std(hv_v, 'omitnan');
    igd_mean(ci)  = mean(igd_v, 'omitnan');
    igd_std(ci)   = std(igd_v, 'omitnan');
    feas_mean(ci) = mean(feas_v, 'omitnan');
    E_mean(ci)    = mean(e_v, 'omitnan');
    runtime_m(ci) = mean(rt_v, 'omitnan');
    if ~isnan(E_baseline)
        ensav_m(ci) = mean(100 * (E_baseline - e_v) / E_baseline, 'omitnan');
    end
end

n_runs_each = max(run_count);

slack_pcts = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100];
T_thresholds = T_min_est * (1 + slack_pcts/100);
n_thr = numel(T_thresholds);
E_thr_mean = nan(n_cfg, n_thr);
E_thr_std  = nan(n_cfg, n_thr);
for ci = 1:n_cfg
    cfg_runs = results(strcmp({results.config_id}, cfg_ids{ci}));
    for ti = 1:n_thr
        E_run = nan(numel(cfg_runs),1);
        for ri = 1:numel(cfg_runs)
            F_here = local_valid_front(cfg_runs(ri));
            if isempty(F_here), continue; end
            feas = F_here(F_here(:,1) <= T_thresholds(ti), :);
            if ~isempty(feas)
                E_run(ri) = min(feas(:,2));
            end
        end
        E_thr_mean(ci,ti) = mean(E_run, 'omitnan');
        E_thr_std(ci,ti)  = std(E_run, 'omitnan');
    end
end

fprintf('Segment  : %s\n', segment_name);
fprintf('MAT file : %s\n', mat_file);
fprintf('T_min    : %.2f s\n', T_min_est);
fprintf('T_sched  : %.2f s\n', T_sched);
fprintf('Runs     : %d per config\n', n_runs_each);
if ~isnan(E_baseline)
    fprintf('E_base   : %.4f kWh\n', E_baseline);
end
fprintf('\n%-8s  %8s  %8s  %8s  %8s  %8s\n', 'Config', 'HV', 'HV_std', 'IGD', 'Feas', 'E_best');
for ci = 1:n_cfg
    fprintf('%-8s  %8.4f  %8.4f  %8.4f  %8.3f  %8.4f\n', ...
        cfg_ids{ci}, hv_mean(ci), hv_std(ci), igd_mean(ci), feas_mean(ci), E_mean(ci));
end
fprintf('\n');

figure(5101); clf;
set(gcf, 'Name', sprintf('E5 HV — %s', segment_name), 'NumberTitle', 'off', ...
    'Position', [50 580 560 430]);
b = bar(hv_mean, 0.65, 'FaceColor', 'flat');
hold on;
errorbar(1:n_cfg, hv_mean, hv_std, 'k.', 'LineWidth', 1.8, 'CapSize', 10);
for ci = 1:n_cfg
    b.CData(ci,:) = CLR.(cfg_ids{ci});
    text(ci, hv_mean(ci) + hv_std(ci) + max(hv_std,[],'omitnan')*0.35 + 1e-3, ...
        sprintf('%.4f', hv_mean(ci)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'FontWeight', 'bold');
end
set(gca, 'XTick', 1:n_cfg, 'XTickLabel', cfg_ids, 'FontSize', 10);
ylabel('Hypervolume (HV)');
title(sprintf('E5: Hypervolume Comparison — %s  (%d runs/config)', segment_name, n_runs_each));
grid on; box on;

figure(5102); clf;
set(gcf, 'Name', sprintf('E5 IGD — %s', segment_name), 'NumberTitle', 'off', ...
    'Position', [640 580 560 430]);
b = bar(igd_mean, 0.65, 'FaceColor', 'flat');
hold on;
errorbar(1:n_cfg, igd_mean, igd_std, 'k.', 'LineWidth', 1.8, 'CapSize', 10);
for ci = 1:n_cfg
    b.CData(ci,:) = CLR.(cfg_ids{ci});
    text(ci, igd_mean(ci) + igd_std(ci) + max(igd_std,[],'omitnan')*0.35 + 1e-3, ...
        sprintf('%.4f', igd_mean(ci)), 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'FontWeight', 'bold');
end
set(gca, 'XTick', 1:n_cfg, 'XTickLabel', cfg_ids, 'FontSize', 10);
ylabel('IGD');
title(sprintf('E5: IGD Comparison — %s  (%d runs/config)', segment_name, n_runs_each));
grid on; box on;

figure(5103); clf;
set(gcf, 'Name', sprintf('E5 Pareto Front — %s', segment_name), 'NumberTitle', 'off', ...
    'Position', [50 80 880 560]);
axes('Position', [0.10 0.28 0.86 0.62]);
hold on;
h_scatter = gobjects(n_cfg, 1);
for ci = 1:n_cfg
    F_union = [];
    cfg_runs = results(strcmp({results.config_id}, cfg_ids{ci}));
    for ri = 1:numel(cfg_runs)
        F_here = local_valid_front(cfg_runs(ri));
        if ~isempty(F_here)
            F_union = [F_union; F_here]; %#ok<AGROW>
        end
    end
    if isempty(F_union), continue; end
    F_nd = nd_filter_local(F_union);
    h_scatter(ci) = scatter(F_nd(:,1), F_nd(:,2), 32, CLR.(cfg_ids{ci}), 'filled', ...
        'Marker', MRK.(cfg_ids{ci}), 'MarkerFaceAlpha', 0.72, ...
        'DisplayName', sprintf('%s  [T %.1f-%.1f s]', LBL.(cfg_ids{ci}), min(F_nd(:,1)), max(F_nd(:,1))));
end
if ~isempty(T_sched) && ~isnan(T_sched)
    xline(T_sched, '-.r', 'LineWidth', 1.5, 'Label', sprintf('T_{sched} = %.0f s', T_sched), ...
        'LabelVerticalAlignment', 'top', 'FontSize', 9);
end
xlabel('Travel time T (s)');
ylabel('Energy consumption E (kWh)');
title(sprintf('E5: Pareto Front Overlay — %s', segment_name));
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

btn_w = 0.14;
btn_h = 0.06;
btn_gap = 0.015;
total_w = n_cfg * btn_w + (n_cfg - 1) * btn_gap;
x_start = (1 - total_w) / 2;
for ci = 1:n_cfg
    clr = CLR.(cfg_ids{ci});
    x_btn = x_start + (ci - 1) * (btn_w + btn_gap);
    uicontrol('Style', 'togglebutton', 'Units', 'normalized', ...
        'Position', [x_btn, 0.03, btn_w, btn_h], 'Value', 1, ...
        'String', sprintf('%s ON/OFF', cfg_ids{ci}), ...
        'BackgroundColor', clr, ...
        'ForegroundColor', ternary_str(mean(clr) < 0.65, 'w', 'k'), ...
        'FontSize', 9, 'FontWeight', 'bold', ...
        'Callback', @(src,~) set(h_scatter(ci), 'Visible', ternary_str(src.Value == 1, 'on', 'off')));
end

figure(5104); clf;
set(gcf, 'Name', sprintf('E5 Best Energy vs Time Threshold — %s', segment_name), ...
    'NumberTitle', 'off', 'Position', [970 80 860 560]);
ax4 = axes('Position', [0.10 0.34 0.86 0.54]);
hold on;
h_line = gobjects(n_cfg, 1);
h_fill = gobjects(n_cfg, 1);
x_labels = cell(1, n_thr);
for ti = 1:n_thr
    x_labels{ti} = sprintf('%.1f s  (+%.0f%%)', T_thresholds(ti), slack_pcts(ti));
end
for ci = 1:n_cfg
    id = cfg_ids{ci};
    % D_LLM is the proposed: thicker line + larger marker, like D in E1 Fig 4
    is_hero = strcmp(id, 'D_LLM');
    lw = 1.5 + is_hero * 1.5;
    ms = 5   + is_hero * 3;
    h_line(ci) = plot(1:n_thr, E_thr_mean(ci,:), ['-' MRK.(id)], ...
        'Color', CLR.(id), 'LineWidth', lw, 'MarkerSize', ms, ...
        'MarkerFaceColor', CLR.(id), 'DisplayName', LBL.(id));
    x_fill = [1:n_thr, fliplr(1:n_thr)];
    y_fill = [E_thr_mean(ci,:) + E_thr_std(ci,:), fliplr(E_thr_mean(ci,:) - E_thr_std(ci,:))];
    valid = ~isnan(y_fill);
    if sum(valid) > 2
        h_fill(ci) = fill(x_fill(valid), y_fill(valid), CLR.(id), ...
            'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
if ~isempty(T_sched) && ~isnan(T_sched)
    x_sched = local_interp_tick(T_thresholds, T_sched);
    xline(x_sched, '-.', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2, ...
        'Label', sprintf('T_{sched} = %.0f s', T_sched), 'HandleVisibility', 'off');
    plot(NaN, NaN, '-.', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('T_{sched} = %.0f s  (+%.0f%% above T_{min})', ...
        T_sched, (T_sched - T_min_est)/T_min_est*100));
end
set(ax4, 'XTick', 1:n_thr, 'XTickLabel', x_labels, 'FontSize', 9, ...
    'TickLabelInterpreter', 'tex');
xtickangle(35);
xlabel(sprintf('Time threshold  (T_{min} = %.2f s — minimum possible travel time)', T_min_est), ...
    'FontSize', 11);
ylabel('Best feasible energy E_{best} (kWh)', 'FontSize', 12);
title(sprintf('E5: Best Energy at Various Time Thresholds — %s  (%d runs/config)', ...
    segment_name, n_runs_each), 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

if ~isempty(T_sched) && ~isnan(T_sched)
    cap_sched = sprintf(' T_{sched} = %.0f s (+%d%% above T_{min}).', ...
        T_sched, round((T_sched - T_min_est)/T_min_est*100));
else
    cap_sched = '';
end
annotation('textbox', [0.05 0.14 0.90 0.07], ...
    'String', ['E_{best} = lowest energy among Pareto solutions where T \leq threshold. ' ...
    sprintf('T_{min} = %.2f s is the fastest possible travel time. ', T_min_est) ...
    'Each threshold is shown as actual time (s) and its % above T_{min}.' ...
    cap_sched ' Shaded bands = \pm1 std across runs. D_LLM (bold) = proposed + LLM advisor.'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 8.5, 'FontAngle', 'italic', 'VerticalAlignment', 'middle');

btn_w = 0.14;
btn_h = 0.06;
btn_gap = 0.015;
total_w = n_cfg * btn_w + (n_cfg - 1) * btn_gap;
x_start = (1 - total_w) / 2;
for ci = 1:n_cfg
    clr = CLR.(cfg_ids{ci});
    hi = h_fill(ci);
    x_btn = x_start + (ci - 1) * (btn_w + btn_gap);
    uicontrol('Style', 'togglebutton', 'Units', 'normalized', ...
        'Position', [x_btn, 0.03, btn_w, btn_h], 'Value', 1, ...
        'String', sprintf('%s ON/OFF', cfg_ids{ci}), ...
        'BackgroundColor', clr, ...
        'ForegroundColor', ternary_str(mean(clr) < 0.65, 'w', 'k'), ...
        'FontSize', 9, 'FontWeight', 'bold', ...
        'Callback', @(src,~) toggle_line(src, h_line(ci), hi));
end

%% =================================================================
%% FIGURE 6 — Energy Comparison Table: D_LLM vs all others
%% =================================================================
figure(5106); clf;
set(gcf, 'Name', sprintf('E5 Energy Comparison Table — %s', segment_name), ...
    'NumberTitle', 'off', 'Position', [50 80 1200 400]);

sel_slacks = [10, 20, 30, 40, 50, 70];
sel_T      = T_min_est * (1 + sel_slacks/100);
if ~isempty(T_sched) && ~isnan(T_sched)
    sel_slacks(end+1) = round((T_sched - T_min_est)/T_min_est*100);
    sel_T(end+1)      = T_sched;
end
n_sel = numel(sel_T);

% Compute E_best per config per selected threshold
E_tbl6 = nan(n_cfg, n_sel);
for ci = 1:n_cfg
    cfg_runs = results(strcmp({results.config_id}, cfg_ids{ci}));
    for ti = 1:n_sel
        E_run = nan(numel(cfg_runs),1);
        for ri = 1:numel(cfg_runs)
            F_here = local_valid_front(cfg_runs(ri));
            if isempty(F_here), continue; end
            fp = F_here(F_here(:,1) <= sel_T(ti), :);
            if ~isempty(fp), E_run(ri) = min(fp(:,2)); end
        end
        E_tbl6(ci,ti) = mean(E_run, 'omitnan');
    end
end

% D_LLM index
idx_DLLM = find(strcmp(cfg_ids, 'D_LLM'));
compare_ids = cfg_ids(~strcmp(cfg_ids, 'D_LLM'));

col_names = [{'T threshold (s)', 'Slack'}, ...
    cellfun(@(x) sprintf('E_%s (kWh)', x), cfg_ids, 'UniformOutput', false), ...
    cellfun(@(x) sprintf('D_{LLM} vs %s (kWh)', x), compare_ids, 'UniformOutput', false), ...
    cellfun(@(x) sprintf('D_{LLM} vs %s (%%)', x), compare_ids, 'UniformOutput', false)];

tbl_data6 = cell(n_sel, numel(col_names));
for ti = 1:n_sel
    is_sched = ~isempty(T_sched) && ~isnan(T_sched) && (ti == n_sel);
    T_label  = sprintf('%.1f%s', sel_T(ti), ternary_str(is_sched, ' ★', ''));
    E_DLLM   = E_tbl6(idx_DLLM, ti);

    row = {T_label, sprintf('+%.0f%%', sel_slacks(ti))};
    for ci = 1:n_cfg
        row{end+1} = sprintf('%.4f', E_tbl6(ci,ti)); %#ok<AGROW>
    end
    for ci = 1:n_cfg
        if strcmp(cfg_ids{ci}, 'D_LLM'), continue; end
        E_other = E_tbl6(ci, ti);
        row{end+1} = sprintf('%+.4f', E_other - E_DLLM); %#ok<AGROW>
    end
    for ci = 1:n_cfg
        if strcmp(cfg_ids{ci}, 'D_LLM'), continue; end
        E_other = E_tbl6(ci, ti);
        row{end+1} = sprintf('%+.2f%%', (E_other - E_DLLM) / E_other * 100); %#ok<AGROW>
    end
    tbl_data6(ti,:) = row;
end

uitable('Data', tbl_data6, 'ColumnName', col_names, ...
    'Units', 'normalized', 'Position', [0.01 0.15 0.98 0.75], ...
    'FontSize', 9, 'RowName', {});

annotation('textbox', [0.01 0.88 0.98 0.10], ...
    'String', sprintf(['E5 Energy Comparison: D_LLM vs All Configs  |  Segment %s  |  ' ...
    'T_{min} = %.2f s  |  %d runs/config'], segment_name, T_min_est, n_runs_each), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'middle');
annotation('textbox', [0.01 0.01 0.98 0.10], ...
    'String', ['Positive value = D_LLM uses LESS energy (D_LLM is better). ' ...
    'E_{best} = lowest energy among Pareto solutions with T \leq T_{threshold}. ' ...
    '★ = actual scheduled time (T_{sched}).'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 8.5, 'FontAngle', 'italic');

% Print table to console
fprintf('\n=== Energy Comparison: D_LLM vs Others | %s ===\n', segment_name);
fprintf('%-12s  %-6s', 'T_thr(s)', 'Slack');
for ci = 1:n_cfg, fprintf('  %10s', sprintf('E_%s', cfg_ids{ci})); end
for ci = 1:n_cfg
    if ~strcmp(cfg_ids{ci},'D_LLM'), fprintf('  %+10s', sprintf('vs%s(kWh)',cfg_ids{ci})); end
end
fprintf('\n%s\n', repmat('-',1,130));
for ti = 1:n_sel
    mk = ''; if ~isempty(T_sched) && ~isnan(T_sched) && ti==n_sel, mk=' ★'; end
    fprintf('%-12s  %-6s', sprintf('%.1f%s',sel_T(ti),mk), sprintf('+%.0f%%',sel_slacks(ti)));
    for ci = 1:n_cfg, fprintf('  %10.4f', E_tbl6(ci,ti)); end
    for ci = 1:n_cfg
        if ~strcmp(cfg_ids{ci},'D_LLM')
            fprintf('  %+10.4f', E_tbl6(ci,ti) - E_tbl6(idx_DLLM,ti));
        end
    end
    fprintf('\n');
end
fprintf('\n');

figure(5105); clf;
set(gcf, 'Name', sprintf('E5 Summary Table — %s', segment_name), 'NumberTitle', 'off', ...
    'Position', [120 120 1120 380]);
tbl = cell(n_cfg, 8);
for ci = 1:n_cfg
    tbl(ci,:) = {cfg_ids{ci}, sprintf('%.4f', hv_mean(ci)), sprintf('%.4f', hv_std(ci)), ...
        sprintf('%.4f', igd_mean(ci)), sprintf('%.3f', feas_mean(ci)), ...
        sprintf('%.4f', E_mean(ci)), sprintf('%.2f', ensav_m(ci)), sprintf('%.1f', runtime_m(ci))};
end
uitable('Data', tbl, ...
    'ColumnName', {'Config','HV mean','HV std','IGD mean','Feas mean','E_best mean','EnSav %','Runtime mean (s)'}, ...
    'Units', 'normalized', 'Position', [0.02 0.14 0.96 0.78], 'FontSize', 10, 'RowName', {});
annotation('textbox', [0.02 0.90 0.96 0.08], ...
    'String', sprintf('E5 Summary Reloaded from MAT — %s | T_{min}=%.2f s | T_{sched}=%.2f s | Runs=%d', ...
    segment_name, T_min_est, T_sched, n_runs_each), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
annotation('textbox', [0.02 0.02 0.96 0.08], ...
    'String', 'If EnSav is NaN, the MAT file did not include baseline metadata and baseline.mat was not found.', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8.5, 'FontAngle', 'italic');

fprintf('Done — 6 figures displayed from MAT data.\n');
end

function out = local_get_field(s, name, default)
    if isstruct(s) && isfield(s, name)
        out = s.(name);
    else
        out = default;
    end
end

function txt = local_pick_text(default_txt, varargin)
    txt = default_txt;
    for i = 1:numel(varargin)
        cand = varargin{i};
        if ischar(cand) || (isstring(cand) && isscalar(cand))
            cand = char(cand);
            if ~isempty(strtrim(cand))
                txt = cand;
                return;
            end
        end
    end
end

function F_here = local_valid_front(r)
    F_here = zeros(0,2);
    if isfield(r, 'F') && ~isempty(r.F)
        F_here = double(r.F);
        return;
    end
    if isfield(r, 'AllF') && ~isempty(r.AllF)
        F_all = double(r.AllF);
        valid = F_all(:,1) < 1e5 & F_all(:,2) < 1e5;
        F_here = nd_filter_local(F_all(valid,:));
    end
end

function x_tick = local_interp_tick(T_thresholds, T_sched)
    if T_sched <= T_thresholds(1)
        x_tick = 1;
    elseif T_sched >= T_thresholds(end)
        x_tick = numel(T_thresholds);
    else
        idx_lo = find(T_thresholds <= T_sched, 1, 'last');
        idx_hi = idx_lo + 1;
        frac = (T_sched - T_thresholds(idx_lo)) / (T_thresholds(idx_hi) - T_thresholds(idx_lo));
        x_tick = idx_lo + frac;
    end
end

function s = ternary_str(cond, a, b)
    if cond
        s = a;
    else
        s = b;
    end
end

function toggle_line(btn, h_line, h_fill)
    vis = ternary_str(btn.Value == 1, 'on', 'off');
    set(h_line, 'Visible', vis);
    if isgraphics(h_fill)
        set(h_fill, 'Visible', vis);
    end
end

function F_nd = nd_filter_local(F)
    if isempty(F)
        F_nd = zeros(0,2);
        return;
    end
    n = size(F,1);
    keep = true(n,1);
    for i = 1:n
        if ~keep(i), continue; end
        for j = 1:n
            if i == j || ~keep(j), continue; end
            if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                keep(i) = false;
                break;
            end
        end
    end
    F_nd = F(keep,:);
end
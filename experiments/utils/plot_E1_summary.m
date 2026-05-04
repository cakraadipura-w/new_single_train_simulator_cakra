function plot_E1_summary(seg_name, T_sched)
%% PLOT_E1_SUMMARY  Visualize E1 ablation results from saved experiment files.
%
% Generates 4 figures (all labels in English):
%   Fig 1 — Hypervolume (HV) bar chart per config
%   Fig 2 — IGD bar chart per config
%   Fig 3 — Pareto front overlay: A, B, C, D (bold colors)
%   Fig 4 — Best energy at multiple T-thresholds (various slack levels)
%
% Usage:
%   plot_E1_summary                    % default: IS02
%   plot_E1_summary('IS04')
%   plot_E1_summary('IS02', 170)       % provide scheduled time T_sched (s)

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

if nargin < 1 || isempty(seg_name), seg_name = 'IS07'; end
if nargin < 2, T_sched = []; end   % optional: actual timetable T (s)

% Auto-fill T_sched from known segment schedule times (from main_experiments.m)
T_SCHED_LUT = struct('IS01',130,'IS02',170,'IS03',185,'IS04',180, ...
                     'IS05',185,'IS06',220,'IS07',210,'IS08',330);
if isempty(T_sched) && isfield(T_SCHED_LUT, seg_name)
    T_sched = T_SCHED_LUT.(seg_name);
end

%% --- Load data files ---
res_dir  = fullfile(pwd, 'experiment_results', seg_name);
csv_file = fullfile(res_dir, 'E1_results.csv');
mat_file = fullfile(res_dir, 'E1_all_configs.mat');

assert(exist(csv_file,'file')==2, 'File not found: %s', csv_file);
assert(exist(mat_file,'file')==2, 'File not found: %s', mat_file);

T_csv   = readtable(csv_file);
data    = load(mat_file);
results = data.results;

%% --- Config colors and styles (distinct, easy to tell apart) ---
% A = gray (baseline), B = blue, C = orange, D = bright green (proposed)
CLR = struct('A',[0.55 0.55 0.55], 'B',[0.20 0.45 0.85], ...
             'C',[0.95 0.50 0.05], 'D',[0.08 0.78 0.18]);
MRK = struct('A','o',  'B','s',  'C','^',  'D','d');   % marker shapes
LBL = struct( ...
    'A','A: Baseline (Orig. CC\_CR + Orig. NSGA-II)', ...
    'B','B: Solver only (Orig. CC\_CR + BRL-SDE)', ...
    'C','C: Strategy only (Impr. CC\_CR + Orig. NSGA-II)', ...
    'D','D: Proposed (Impr. CC\_CR + BRL-SDE)');

cfg_ids = {'A','B','C','D'};
n_cfg   = numel(cfg_ids);

%% --- Estimate T_min from Pareto data ---
all_T = [];
for i = 1:numel(results)
    if ~isempty(results(i).F)
        all_T = [all_T; results(i).F(:,1)]; %#ok<AGROW>
    end
end
T_min_est = min(all_T);
n_runs_each = sum(strcmp(T_csv.config, 'A'));

fprintf('Segment  : %s\n', seg_name);
fprintf('T_min    : %.2f s  (estimated from Pareto data)\n', T_min_est);
if ~isempty(T_sched)
    fprintf('T_sched  : %.1f s  (slack %.1f%% above T_min)\n', ...
        T_sched, (T_sched - T_min_est)/T_min_est*100);
end
fprintf('Runs     : %d per config\n\n', n_runs_each);

%% --- Collect HV / IGD statistics from CSV ---
hv_mean  = zeros(n_cfg,1);  hv_std  = zeros(n_cfg,1);
igd_mean = zeros(n_cfg,1);  igd_std = zeros(n_cfg,1);

for ci = 1:n_cfg
    mask = strcmp(T_csv.config, cfg_ids{ci});
    hv_v  = T_csv.hv(mask);
    igd_v = T_csv.igd(mask);
    hv_mean(ci)  = mean(hv_v,  'omitnan');
    hv_std(ci)   = std(hv_v,   'omitnan');
    igd_mean(ci) = mean(igd_v, 'omitnan');
    igd_std(ci)  = std(igd_v,  'omitnan');
end

%% --- Compute E_best at multiple T-thresholds ---
% Range: 5% to 100% above T_min — covers before AND after T_sched
slack_pcts   = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100];
T_thresholds = T_min_est * (1 + slack_pcts/100);
% T_sched is NOT added as a data point — it is marked as a vertical line
% at its interpolated x-axis position (see Figure 4)

n_thr = numel(T_thresholds);
E_mean = nan(n_cfg, n_thr);   % mean E_best per config per threshold
E_std  = nan(n_cfg, n_thr);   % std  E_best

for ci = 1:n_cfg
    runs = find(strcmp({results.config_id}, cfg_ids{ci}));
    for ti = 1:n_thr
        T_thr = T_thresholds(ti);
        E_run = nan(numel(runs),1);
        for ri = 1:numel(runs)
            F = results(runs(ri)).F;
            if isempty(F), continue; end
            fp = F(F(:,1) <= T_thr, :);
            if ~isempty(fp), E_run(ri) = min(fp(:,2)); end
        end
        E_mean(ci,ti) = mean(E_run, 'omitnan');
        E_std(ci,ti)  = std(E_run,  'omitnan');
    end
end

%% --- Print summary table ---
fprintf('%-8s  %8s  %8s  %8s  %8s\n','Config','HV','HV_std','IGD','IGD_std');
for ci = 1:n_cfg
    fprintf('%-8s  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
        cfg_ids{ci}, hv_mean(ci), hv_std(ci), igd_mean(ci), igd_std(ci));
end
fprintf('\n');

%% =================================================================
%% FIGURE 1 — Hypervolume (HV) Bar Chart
%% =================================================================
figure(5001); clf;
set(gcf,'Name',sprintf('E1 HV — %s',seg_name),'NumberTitle','off', ...
    'Position',[50 580 520 430]);

ax1 = axes('Position',[0.13 0.22 0.82 0.65]);
b = bar(hv_mean, 0.6, 'FaceColor','flat', 'Parent', ax1);
for ci = 1:n_cfg, b.CData(ci,:) = CLR.(cfg_ids{ci}); end
hold on;
errorbar(1:n_cfg, hv_mean, hv_std, 'k.', 'LineWidth',2, 'CapSize',10);

% Value labels on top of each bar
y_offset = range(hv_mean)*0.04 + max(hv_std)*0.5;
for ci = 1:n_cfg
    text(ci, hv_mean(ci)+hv_std(ci)+y_offset, ...
        sprintf('%.4f', hv_mean(ci)), ...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
end

set(ax1,'XTick',1:n_cfg,'XTickLabel',cfg_ids,'FontSize',12);
ylabel('Hypervolume (HV)  ↑ higher is better','FontSize',12);
title(sprintf('E1: Hypervolume Comparison — %s  (%d runs/config)', ...
    seg_name, n_runs_each),'FontSize',12);
ylim([min(hv_mean)-range(hv_mean)*0.5, max(hv_mean)+range(hv_mean)*0.8]);
grid on; box on;

annotation('textbox',[0.05 0.01 0.90 0.10], ...
    'String',['Higher HV = better overall solution quality. ' ...
    'D (green) indicates the proposed method covers a larger portion of the ' ...
    'Pareto-optimal objective space than the baselines.'], ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.5,'FontAngle','italic','VerticalAlignment','middle');


%% =================================================================
%% FIGURE 2 — IGD Bar Chart
%% =================================================================
figure(5002); clf;
set(gcf,'Name',sprintf('E1 IGD — %s',seg_name),'NumberTitle','off', ...
    'Position',[590 580 520 430]);

ax2 = axes('Position',[0.13 0.22 0.82 0.65]);
b = bar(igd_mean, 0.6, 'FaceColor','flat', 'Parent', ax2);
for ci = 1:n_cfg, b.CData(ci,:) = CLR.(cfg_ids{ci}); end
hold on;
errorbar(1:n_cfg, igd_mean, igd_std, 'k.', 'LineWidth',2, 'CapSize',10);

y_offset = range(igd_mean)*0.04 + max(igd_std)*0.5;
for ci = 1:n_cfg
    text(ci, igd_mean(ci)+igd_std(ci)+y_offset, ...
        sprintf('%.4f', igd_mean(ci)), ...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
end

set(ax2,'XTick',1:n_cfg,'XTickLabel',cfg_ids,'FontSize',12);
ylabel('IGD (lower = better)','FontSize',12);
title(sprintf('E1: IGD Comparison — %s  (%d runs/config)', ...
    seg_name, n_runs_each),'FontSize',12);
ylim([0, max(igd_mean)+range(igd_mean)*0.8]);
grid on; box on;

annotation('textbox',[0.05 0.01 0.90 0.10], ...
    'String',['Lower IGD = Pareto front is closer to the global reference front. ' ...
    'Config D achieves competitive IGD, confirming improved convergence ' ...
    'when using both the enhanced driving strategy and adaptive solver.'], ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.5,'FontAngle','italic','VerticalAlignment','middle');


%% =================================================================
%% FIGURE 3 — Pareto Front Overlay
%% =================================================================
figure(5003); clf;
set(gcf,'Name',sprintf('E1 Pareto Front — %s',seg_name),'NumberTitle','off', ...
    'Position',[50 80 820 560]);

% Axes — leave space at bottom for toggle buttons
axes('Position',[0.10 0.28 0.86 0.62]);
hold on;

% Dot size: D larger so it stands out
SZ = struct('A',28, 'B',28, 'C',28, 'D',50);

% Print T range per config to console
fprintf('\n--- Pareto Front T range per Config (all runs combined) ---\n');
fprintf('%-8s  %10s  %10s  %10s  %8s\n', ...
    'Config', 'T_min(s)', 'T_max(s)', 'E_min(kWh)', 'N_pts');

h_scatter = gobjects(n_cfg, 1);   % store scatter handles for toggle buttons

for ci = 1:n_cfg
    id   = cfg_ids{ci};
    runs = find(strcmp({results.config_id}, id));
    F_union = [];
    for ri = 1:numel(runs)
        F = results(runs(ri)).F;
        if ~isempty(F), F_union = [F_union; F]; end %#ok<AGROW>
    end
    if isempty(F_union), continue; end
    F_nd = nd_filter_local(F_union);

    T_lo = min(F_nd(:,1));
    T_hi = max(F_nd(:,1));
    E_lo = min(F_nd(:,2));

    fprintf('%-8s  %10.2f  %10.2f  %10.4f  %8d\n', ...
        id, T_lo, T_hi, E_lo, size(F_nd,1));

    h_scatter(ci) = scatter(F_nd(:,1), F_nd(:,2), SZ.(id), CLR.(id), 'filled', ...
        'MarkerFaceAlpha', 0.75, ...
        'DisplayName', sprintf('%s  [T: %.1f–%.1f s]', LBL.(id), T_lo, T_hi));
end
fprintf('\n');

% Mark T_sched if provided
if ~isempty(T_sched)
    xline(T_sched, '-.r', 'LineWidth',1.5, ...
        'Label', sprintf('T_{sched} = %.0fs', T_sched), ...
        'LabelVerticalAlignment','top','FontSize',9);
end

xlabel('Travel time T (s)','FontSize',12);
ylabel('Energy consumption E (kWh)','FontSize',12);
title(sprintf('E1: Pareto Front Overlay — %s  (all runs combined)', seg_name),'FontSize',12);
legend('Location','northeast','FontSize',9);
grid on; box on;
xlim([T_min_est*0.97, 400]);

% ---- Toggle buttons (one per config) ----
% Click to show/hide each config's Pareto points
btn_w   = 0.20;
btn_h   = 0.06;
btn_gap = 0.025;
total_w = n_cfg * btn_w + (n_cfg-1) * btn_gap;
x_start = (1 - total_w) / 2;   % center the buttons

for ci = 1:n_cfg
    id   = cfg_ids{ci};
    clr  = CLR.(id);
    % Use white text on dark backgrounds, black on light
    txt_clr = ternary_str(mean(clr) < 0.65, 'w', 'k');
    x_btn = x_start + (ci-1) * (btn_w + btn_gap);

    uicontrol('Style','togglebutton', ...
        'Units','normalized', ...
        'Position',[x_btn, 0.03, btn_w, btn_h], ...
        'String', sprintf('Config %s  ON/OFF', id), ...
        'Value', 1, ...
        'BackgroundColor', clr, ...
        'ForegroundColor', txt_clr, ...
        'FontSize', 10, 'FontWeight','bold', ...
        'Callback', @(src,~) set(h_scatter(ci), ...
            'Visible', ternary_str(src.Value==1,'on','off')));
end

annotation('textbox',[0.05 0.10 0.90 0.06], ...
    'String',['Click the buttons below to show/hide each config. ' ...
    'Points in the lower-left are better (fast & energy-efficient). ' ...
    'Config D (green, proposed method) should dominate near T_{sched}.'], ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.5,'FontAngle','italic','VerticalAlignment','middle');


%% =================================================================
%% FIGURE 4 — Best Energy at Multiple T-Thresholds
%% =================================================================
figure(5004); clf;
set(gcf,'Name',sprintf('E1 Best Energy vs Time Threshold — %s',seg_name), ...
    'NumberTitle','off','Position',[590 80 820 560]);

ax4 = axes('Position',[0.10 0.34 0.86 0.54]);
hold on;

% X-axis labels: actual time (s) + slack % above T_min
x_labels = cell(1, n_thr);
for ti = 1:n_thr
    x_labels{ti} = sprintf('%.1f s  (slack %.0f%%)', ...
        T_thresholds(ti), slack_pcts(ti));
end

h_line = gobjects(n_cfg, 1);   % store line handles for toggle buttons
h_fill = gobjects(n_cfg, 1);   % store fill handles

for ci = 1:n_cfg
    id = cfg_ids{ci};
    lw = 1.5 + (strcmp(id,'D')*1.5);
    ms = 5   + (strcmp(id,'D')*3);

    h_line(ci) = plot(1:n_thr, E_mean(ci,:), ['-' MRK.(id)], ...
        'Color', CLR.(id), 'LineWidth', lw, ...
        'MarkerSize', ms, 'MarkerFaceColor', CLR.(id), ...
        'DisplayName', LBL.(id));

    x_fill = [1:n_thr, fliplr(1:n_thr)];
    y_fill = [E_mean(ci,:)+E_std(ci,:), fliplr(E_mean(ci,:)-E_std(ci,:))];
    valid  = ~isnan(y_fill);
    if sum(valid) > 2
        h_fill(ci) = fill(x_fill(valid), y_fill(valid), CLR.(id), ...
            'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
    end
end

% --- Mark T_sched at its interpolated x position (between ticks) ---
if ~isempty(T_sched)
    sched_slack = (T_sched - T_min_est) / T_min_est * 100;
    if T_sched <= T_thresholds(1)
        x_sched = 1;
    elseif T_sched >= T_thresholds(end)
        x_sched = n_thr;
    else
        idx_lo  = find(T_thresholds <= T_sched, 1, 'last');
        idx_hi  = idx_lo + 1;
        frac    = (T_sched - T_thresholds(idx_lo)) / ...
                  (T_thresholds(idx_hi) - T_thresholds(idx_lo));
        x_sched = idx_lo + frac;
    end
    xline(x_sched, '-.', 'Color',[0.85 0.10 0.10], 'LineWidth',2.2, ...
        'Label',sprintf('T_{sched} = %.0f s (+%.0f%%)', T_sched, sched_slack), ...
        'LabelVerticalAlignment','top','FontSize',10, 'HandleVisibility','off');
end

% --- Dummy legend entry for T_sched line ---
if ~isempty(T_sched)
    plot(NaN, NaN, '-.', 'Color',[0.85 0.10 0.10], 'LineWidth',2.2, ...
        'DisplayName', sprintf('Scheduled time  (T_{sched} = %.0f s, +%.0f%% above T_{min})', ...
        T_sched, (T_sched-T_min_est)/T_min_est*100));
end

set(ax4, 'XTick',1:n_thr, 'XTickLabel',x_labels, 'FontSize',9, ...
    'TickLabelInterpreter','tex');
xtickangle(35);
xlabel(sprintf('Time threshold  (T_{min} = %.2f s — minimum possible travel time)', T_min_est), ...
    'FontSize',11);
ylabel('Best energy E_{best} (kWh)','FontSize',12);
title(sprintf('E1: Best Energy at Various Time Thresholds — %s  (%d runs/config)', ...
    seg_name, n_runs_each),'FontSize',12);
legend('Location','northeast','FontSize',9);
grid on; box on;

% Build caption with T_min and T_sched info
if ~isempty(T_sched)
    cap_sched = sprintf(' T_{sched} = %.0f s (+%d%% above T_{min}).', ...
        T_sched, round((T_sched-T_min_est)/T_min_est*100));
else
    cap_sched = '';
end
annotation('textbox',[0.05 0.14 0.90 0.07], ...
    'String',['E_{best} = lowest energy among Pareto solutions where T \leq threshold. ' ...
    sprintf('T_{min} = %.2f s is the fastest possible time (full-throttle, no energy saving). ', T_min_est) ...
    'Each threshold is shown as actual time (s) and its % above T_{min}.' ...
    cap_sched ' Shaded bands = ±1 std across runs.'], ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.5,'FontAngle','italic','VerticalAlignment','middle');

% ---- Toggle buttons (one per config) ----
btn_w   = 0.20;
btn_h   = 0.06;
btn_gap = 0.025;
total_w = n_cfg * btn_w + (n_cfg-1) * btn_gap;
x_start = (1 - total_w) / 2;

for ci = 1:n_cfg
    id      = cfg_ids{ci};
    clr     = CLR.(id);
    txt_clr = ternary_str(mean(clr) < 0.65, 'w', 'k');
    x_btn   = x_start + (ci-1) * (btn_w + btn_gap);
    hi      = h_fill(ci);   % fill handle (may be null object)

    uicontrol('Style','togglebutton', ...
        'Units','normalized', ...
        'Position',[x_btn, 0.03, btn_w, btn_h], ...
        'String', sprintf('Config %s  ON/OFF', id), ...
        'Value', 1, ...
        'BackgroundColor', clr, ...
        'ForegroundColor', txt_clr, ...
        'FontSize', 10, 'FontWeight','bold', ...
        'Callback', @(src,~) toggle_line(src, h_line(ci), hi));
end

%% =================================================================
%% FIGURE 5 — Energy Comparison Table: Config D vs A, B, C
%% =================================================================
figure(5005); clf;
set(gcf,'Name',sprintf('E1 Energy Comparison Table — %s',seg_name), ...
    'NumberTitle','off','Position',[50 80 980 400]);

% Selected time thresholds for the table
sel_slacks = [10, 20, 30, 40, 50, 70];
sel_T      = T_min_est * (1 + sel_slacks/100);

% Append T_sched as a separate row if it's not already covered
if ~isempty(T_sched)
    sel_slacks = [sel_slacks, round((T_sched-T_min_est)/T_min_est*100)];
    sel_T      = [sel_T, T_sched];
end
n_sel = numel(sel_T);

% Compute E_best per config per selected threshold
E_tbl = nan(4, n_sel);   % rows: A B C D
for ci = 1:n_cfg
    runs = find(strcmp({results.config_id}, cfg_ids{ci}));
    for ti = 1:n_sel
        E_run = nan(numel(runs),1);
        for ri = 1:numel(runs)
            F = results(runs(ri)).F;
            if isempty(F), continue; end
            fp = F(F(:,1) <= sel_T(ti), :);
            if ~isempty(fp), E_run(ri) = min(fp(:,2)); end
        end
        E_tbl(ci,ti) = mean(E_run,'omitnan');
    end
end

% Build table cell array
% Columns: T_threshold | Slack | E_A | E_B | E_C | E_D | ΔD-A(kWh) | ΔD-A(%) | ΔD-B(kWh) | ΔD-B(%) | ΔD-C(kWh) | ΔD-C(%)
col_names = {'T threshold (s)', 'Slack above T_min', ...
    'Energy Config A (kWh)', 'Energy Config B (kWh)', 'Energy Config C (kWh)', 'Energy Config D (kWh)', ...
    'D vs A (kWh)', 'D vs A (%)', ...
    'D vs B (kWh)', 'D vs B (%)', ...
    'D vs C (kWh)', 'D vs C (%)'};

tbl_data = cell(n_sel, numel(col_names));
for ti = 1:n_sel
    E_A = E_tbl(1,ti);  E_B = E_tbl(2,ti);
    E_C = E_tbl(3,ti);  E_D = E_tbl(4,ti);

    % Positive = D saves more energy vs that config
    dAkwh = E_A - E_D;  dApct = dAkwh / E_A * 100;
    dBkwh = E_B - E_D;  dBpct = dBkwh / E_B * 100;
    dCkwh = E_C - E_D;  dCpct = dCkwh / E_C * 100;

    is_sched = ~isempty(T_sched) && (ti == n_sel);
    T_label  = sprintf('%.1f%s', sel_T(ti), ternary_str(is_sched,' ★',''));

    tbl_data(ti,:) = { ...
        T_label, ...
        sprintf('+%.0f%%', sel_slacks(ti)), ...
        sprintf('%.4f', E_A), ...
        sprintf('%.4f', E_B), ...
        sprintf('%.4f', E_C), ...
        sprintf('%.4f', E_D), ...
        sprintf('%+.4f', dAkwh), sprintf('%+.2f%%', dApct), ...
        sprintf('%+.4f', dBkwh), sprintf('%+.2f%%', dBpct), ...
        sprintf('%+.4f', dCkwh), sprintf('%+.2f%%', dCpct)};
end

% Render uitable
uit = uitable('Data', tbl_data, 'ColumnName', col_names, ...
    'Units','normalized','Position',[0.01 0.15 0.98 0.75], ...
    'FontSize', 10, 'RowName', {});

% Auto column width
col_w = {95, 95, 80, 80, 80, 80, 90, 80, 90, 80, 90, 80};
set(uit, 'ColumnWidth', col_w);

% Title and caption
annotation('textbox',[0.01 0.88 0.98 0.10], ...
    'String', sprintf(['E1 Energy Comparison: Config D (Proposed) vs A, B, C  |  Segment %s  |  ' ...
    'T_{min} = %.2f s  |  %d runs/config'], seg_name, T_min_est, n_runs_each), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',11,'FontWeight','bold','VerticalAlignment','middle');

annotation('textbox',[0.01 0.01 0.98 0.10], ...
    'String',['Positive value = Config D uses LESS energy (D is better). ' ...
    'E_best = lowest energy among Pareto solutions with T \leq T_{threshold}. ' ...
    '★ = actual scheduled time (T_{sched}).'], ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.5,'FontAngle','italic','VerticalAlignment','middle');

% Print same table to console
fprintf('\n=== Energy Comparison Table: Config D vs A, B, C | %s ===\n', seg_name);
fprintf('%-14s  %-8s  %10s  %10s  %10s  %10s  %10s  %8s  %10s  %8s  %10s  %8s\n', ...
    'T_thr(s)','Slack','E_A(kWh)','E_B(kWh)','E_C(kWh)','E_D(kWh)', ...
    'DvsA(kWh)','DvsA(%)', 'DvsB(kWh)','DvsB(%)','DvsC(kWh)','DvsC(%)');
fprintf('%s\n', repmat('-',1,120));
for ti = 1:n_sel
    E_A=E_tbl(1,ti); E_B=E_tbl(2,ti); E_C=E_tbl(3,ti); E_D=E_tbl(4,ti);
    mk = ''; if ~isempty(T_sched) && ti==n_sel, mk=' ★'; end
    fprintf('%-14s  %-8s  %10.4f  %10.4f  %10.4f  %10.4f  %+10.4f  %+7.2f%%  %+10.4f  %+7.2f%%  %+10.4f  %+7.2f%%%s\n', ...
        sprintf('%.1f',sel_T(ti)), sprintf('+%.0f%%',sel_slacks(ti)), ...
        E_A, E_B, E_C, E_D, ...
        E_A-E_D, (E_A-E_D)/E_A*100, ...
        E_B-E_D, (E_B-E_D)/E_B*100, ...
        E_C-E_D, (E_C-E_D)/E_C*100, mk);
end
fprintf('\n');

fprintf('\nDone — 5 figures displayed.\n');
end

function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end

function toggle_line(btn, h_line, h_fill)
    vis = ternary_str(btn.Value == 1, 'on', 'off');
    set(h_line, 'Visible', vis);
    if isvalid(h_fill)
        set(h_fill, 'Visible', vis);
    end
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

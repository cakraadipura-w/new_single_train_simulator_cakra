%% summarize_e5_schedule_energy.m
% Run this file directly from the MATLAB editor.
%
% Summarizes the E5 Pareto fronts under experiment_results/ISxx by picking
% one representative point near each segment schedule time.
% Rule:
%   - minimum energy within targetToleranceS seconds of T_sched
%   - otherwise the closest-time point, tie-broken by lower energy

clearvars -except selected_segments targetToleranceS route_direction

if ~exist('selected_segments', 'var') || isempty(selected_segments)
    selected_segments = 'all';
end
if ~exist('targetToleranceS', 'var') || isempty(targetToleranceS)
    targetToleranceS = 1.0;
end
if ~exist('route_direction', 'var') || isempty(route_direction)
    route_direction = 'up';
end

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(genpath(project_root));

catalog = get_guangzhou_line7_catalog(route_direction);
catalog = resolve_selected_segments(selected_segments, catalog);

summary_rows = repmat(init_summary_row(), 0, 1);
config_rows = repmat(init_config_row(), 0, 1);

for seg_idx = 1:numel(catalog)
    seg = catalog(seg_idx);
    result_file = fullfile(script_dir, seg.name, 'E5_all_configs.mat');

    if ~isfile(result_file)
        summary_rows(end + 1, 1) = build_missing_summary_row(seg, '', 'missing E5 data'); %#ok<AGROW>
        continue;
    end

    data = load(result_file, 'results');
    assert(isfield(data, 'results'), 'Missing results variable in: %s', result_file);
    results = data.results;

    config_ids = unique(string({results.config_id}), 'stable');
    segment_config_rows = repmat(init_config_row(), 0, 1);

    for cfg_idx = 1:numel(config_ids)
        config_id = config_ids(cfg_idx);
        mask = strcmp(string({results.config_id}), config_id);
        config_results = results(mask);
        if isempty(config_results)
            continue;
        end

        front_table = build_union_front_table(config_results);
        representative = choose_representative_point(front_table, seg.T_sched, targetToleranceS);

        row = init_config_row();
        row.segment = string(seg.name);
        row.target_time_s = seg.T_sched;
        row.experiment_id = "E5";
        row.config_id = config_id;
        row.config_desc = string(config_results(1).config_desc);
        row.nsga2_variant = string(config_results(1).nsga2_variant);
        row.use_improved = logical(config_results(1).use_improved);
        row.representative_time_s = representative.time_s;
        row.representative_energy_kWh = representative.energy_kWh;
        row.time_gap_s = representative.time_s - seg.T_sched;
        row.status = string(representative.status);
        row.point_count = height(front_table);
        row.source_file = string(result_file);

        segment_config_rows(end + 1, 1) = row; %#ok<AGROW>
        config_rows(end + 1, 1) = row; %#ok<AGROW>
    end

    if isempty(segment_config_rows)
        summary_rows(end + 1, 1) = build_missing_summary_row(seg, result_file, 'empty compatible E5 results'); %#ok<AGROW>
        continue;
    end

    best_idx = choose_best_config_row(segment_config_rows);
    best_row = segment_config_rows(best_idx);

    summary = init_summary_row();
    summary.segment = string(seg.name);
    summary.target_time_s = seg.T_sched;
    summary.experiment_id = "E5";
    summary.best_config_id = best_row.config_id;
    summary.best_config_desc = best_row.config_desc;
    summary.best_time_s = best_row.representative_time_s;
    summary.best_energy_kWh = best_row.representative_energy_kWh;
    summary.time_gap_s = best_row.time_gap_s;
    summary.status = best_row.status;
    summary.available_configs = numel(segment_config_rows);
    summary.source_file = string(result_file);
    summary_rows(end + 1, 1) = summary; %#ok<AGROW>
end

summary_table = struct2table(summary_rows);
config_table = struct2table(config_rows);

summary_csv = fullfile(script_dir, 'e5_schedule_energy_summary.csv');
config_csv = fullfile(script_dir, 'e5_schedule_energy_by_config.csv');
summary_txt = fullfile(script_dir, 'e5_schedule_energy_summary.txt');

writetable(summary_table, summary_csv);
writetable(config_table, config_csv);
write_summary_text(summary_txt, summary_table, config_table, targetToleranceS);

disp(summary_table(:, {'segment', 'experiment_id', 'target_time_s', 'best_config_id', 'best_time_s', 'best_energy_kWh', 'time_gap_s', 'status'}));
fprintf('Saved summary CSV : %s\n', summary_csv);
fprintf('Saved config CSV  : %s\n', config_csv);
fprintf('Saved summary TXT : %s\n', summary_txt);

function catalog = resolve_selected_segments(selected_segments, catalog)
    if isnumeric(selected_segments)
        run_idx = selected_segments(:)';
    else
        selected_tokens = cellstr(string(selected_segments));
        if numel(selected_tokens) == 1 && strcmpi(selected_tokens{1}, 'all')
            run_idx = 1:numel(catalog);
        else
            run_idx = zeros(1, numel(selected_tokens));
            for idx = 1:numel(selected_tokens)
                seg_idx = find(strcmpi({catalog.name}, selected_tokens{idx}), 1, 'first');
                assert(~isempty(seg_idx), 'Invalid segment token: %s', selected_tokens{idx});
                run_idx(idx) = seg_idx;
            end
        end
    end

    run_idx = unique(run_idx, 'stable');
    assert(~isempty(run_idx), 'Select at least one segment.');
    catalog = catalog(run_idx);
end

function row = init_summary_row()
    row = struct( ...
        'segment', "", ...
        'target_time_s', NaN, ...
        'experiment_id', "", ...
        'best_config_id', "", ...
        'best_config_desc', "", ...
        'best_time_s', NaN, ...
        'best_energy_kWh', NaN, ...
        'time_gap_s', NaN, ...
        'status', "", ...
        'available_configs', 0, ...
        'source_file', "");
end

function row = init_config_row()
    row = struct( ...
        'segment', "", ...
        'target_time_s', NaN, ...
        'experiment_id', "", ...
        'config_id', "", ...
        'config_desc', "", ...
        'nsga2_variant', "", ...
        'use_improved', false, ...
        'representative_time_s', NaN, ...
        'representative_energy_kWh', NaN, ...
        'time_gap_s', NaN, ...
        'status', "", ...
        'point_count', 0, ...
        'source_file', "");
end

function row = build_missing_summary_row(seg, source_file, status_text)
    row = init_summary_row();
    row.segment = string(seg.name);
    row.target_time_s = seg.T_sched;
    row.status = string(status_text);
    row.source_file = string(source_file);
end

function front_table = build_union_front_table(results)
    time_vals = zeros(0, 1);
    energy_vals = zeros(0, 1);

    for idx = 1:numel(results)
        if ~isfield(results(idx), 'F') || isempty(results(idx).F)
            continue;
        end
        F = double(results(idx).F);
        if size(F, 2) < 2
            continue;
        end
        mask = all(isfinite(F(:, 1:2)), 2);
        F = F(mask, 1:2);
        time_vals = [time_vals; F(:, 1)]; %#ok<AGROW>
        energy_vals = [energy_vals; F(:, 2)]; %#ok<AGROW>
    end

    front_table = table(time_vals, energy_vals, 'VariableNames', {'time_s', 'energy_kWh'});
    if isempty(front_table)
        return;
    end

    keep_mask = nondominated_mask(front_table{:, {'time_s', 'energy_kWh'}});
    front_table = front_table(keep_mask, :);
    [~, order] = sortrows([front_table.time_s, front_table.energy_kWh], [1 2]);
    front_table = front_table(order, :);
end

function representative = choose_representative_point(front_table, target_time_s, target_tolerance_s)
    representative = struct('time_s', NaN, 'energy_kWh', NaN, 'status', "missing");
    if isempty(front_table)
        return;
    end

    time_delta = abs(front_table.time_s - target_time_s);
    in_window = time_delta <= target_tolerance_s;
    if any(in_window)
        window_table = front_table(in_window, :);
        rel_idx = best_row_index([window_table.energy_kWh, abs(window_table.time_s - target_time_s)]);
        chosen = window_table(rel_idx, :);
        representative.status = "within target window";
    else
        rel_idx = best_row_index([time_delta, front_table.energy_kWh]);
        chosen = front_table(rel_idx, :);
        representative.status = "closest to target";
    end

    representative.time_s = chosen.time_s(1);
    representative.energy_kWh = chosen.energy_kWh(1);
end

function idx = choose_best_config_row(rows)
    statuses = string({rows.status});
    times = [rows.time_gap_s]';
    energies = [rows.representative_energy_kWh]';
    in_window = statuses == "within target window";

    if any(in_window)
        candidate_idx = find(in_window);
        rel_idx = best_row_index([energies(in_window), abs(times(in_window))]);
        idx = candidate_idx(rel_idx);
    else
        idx = best_row_index([abs(times), energies]);
    end
end

function idx = best_row_index(values)
    [~, order] = sortrows(values);
    idx = order(1);
end

function mask = nondominated_mask(values)
    n = size(values, 1);
    mask = true(n, 1);
    for i = 1:n
        if ~mask(i)
            continue;
        end
        for j = 1:n
            if i == j
                continue;
            end
            if all(values(j, :) <= values(i, :)) && any(values(j, :) < values(i, :))
                mask(i) = false;
                break;
            end
        end
    end
end

function write_summary_text(file_path, summary_table, config_table, target_tolerance_s)
    fid = fopen(file_path, 'w');
    assert(fid ~= -1, 'Failed to open summary text file: %s', file_path);
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'E5 schedule-energy summary\n');
    fprintf(fid, 'Rule: min energy within +/- %.1f s; fallback closest time\n\n', target_tolerance_s);
    fprintf(fid, '%-8s | %-4s | %-10s | %-8s | %-11s | %-13s | %-10s | %-22s\n', ...
        'Segment', 'Exp', 'Target(s)', 'Config', 'Time(s)', 'Energy(kWh)', 'Gap(s)', 'Status');
    fprintf(fid, '%s\n', repmat('-', 1, 103));

    for idx = 1:height(summary_table)
        fprintf(fid, '%-8s | %-4s | %-10.3f | %-8s | %-11.3f | %-13.6f | %-10.3f | %-22s\n', ...
            summary_table.segment(idx), ...
            summary_table.experiment_id(idx), ...
            summary_table.target_time_s(idx), ...
            summary_table.best_config_id(idx), ...
            summary_table.best_time_s(idx), ...
            summary_table.best_energy_kWh(idx), ...
            summary_table.time_gap_s(idx), ...
            summary_table.status(idx));
    end

    fprintf(fid, '\nAll configs\n');
    fprintf(fid, '%-8s | %-4s | %-8s | %-5s | %-11s | %-13s | %-10s | %-22s | %-s\n', ...
        'Segment', 'Exp', 'Config', 'Best?', 'Time(s)', 'Energy(kWh)', 'Gap(s)', 'Status', 'Description');
    fprintf(fid, '%s\n', repmat('-', 1, 130));

    for idx = 1:height(config_table)
        seg = config_table.segment(idx);
        best_cfg = "";
        best_mask = summary_table.segment == seg;
        if any(best_mask)
            best_cfg = summary_table.best_config_id(find(best_mask, 1, 'first'));
        end

        if config_table.config_id(idx) == best_cfg
            best_label = 'yes';
        else
            best_label = '';
        end

        fprintf(fid, '%-8s | %-4s | %-8s | %-5s | %-11.3f | %-13.6f | %-10.3f | %-22s | %s\n', ...
            config_table.segment(idx), ...
            config_table.experiment_id(idx), ...
            config_table.config_id(idx), ...
            best_label, ...
            config_table.representative_time_s(idx), ...
            config_table.representative_energy_kWh(idx), ...
            config_table.time_gap_s(idx), ...
            config_table.status(idx), ...
            config_table.config_desc(idx));
    end
end
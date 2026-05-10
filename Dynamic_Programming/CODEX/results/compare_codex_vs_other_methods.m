%% compare_codex_vs_other_methods.m
% Run this file directly from the MATLAB Editor.
%
% Compares CODEX DP results against:
%   - NSGA-II from main/results
%   - MOPSO from main/results
%   - Convex Optimisation results
%   - MILP results
%   - Improved NSGA-II + improved CC_CR from experiment_results (E1 config D)

clear

if ~exist('selected_segments', 'var') || isempty(selected_segments)
    selected_segments = 'auto';
end

if ~exist('showFigures', 'var') || isempty(showFigures)
    showFigures = true;
end

if ~exist('targetToleranceS', 'var') || isempty(targetToleranceS)
    targetToleranceS = 1.0;
end

if ~exist('preferredCOSolver', 'var') || isempty(preferredCOSolver)
    preferredCOSolver = 'mosek';
end

if ~exist('preferredMILPSolver', 'var') || isempty(preferredMILPSolver)
    preferredMILPSolver = 'mosek';
end

if ~exist('improvedConfigId', 'var') || isempty(improvedConfigId)
    improvedConfigId = 'D';
end

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(fileparts(this_dir)));
addpath(genpath(project_root));

output_root = fullfile(this_dir, 'comparison');
if exist(output_root, 'dir') ~= 7
    mkdir(output_root);
end

summary = compare_codex_vs_other_methods_impl( ...
    project_root, output_root, selected_segments, showFigures, ...
    targetToleranceS, preferredCOSolver, preferredMILPSolver, improvedConfigId);

if ~isempty(summary)
    disp(struct2table(summary));
end

function summary = compare_codex_vs_other_methods_impl(project_root, output_root, selected_segments, showFigures, target_tolerance_s, preferred_co_solver, preferred_milp_solver, improved_config_id)
    codex_bundle = load_codex_bundle(project_root);
    selected_runs = resolve_selected_runs(codex_bundle.segment_runs, selected_segments);

    summary = repmat(struct( ...
        'segment', '', ...
        'segment_csv', '', ...
        'segment_png', '', ...
        'available_methods', 0, ...
        'missing_methods', ''), 0, 1);

    all_rows = table();

    for idx = 1:numel(selected_runs)
        seg_run = selected_runs(idx);
        segment = seg_run.segment;
        segment_name = char(string(segment.name));
        target_time_s = double(segment.T_sched);

        segment_output_dir = fullfile(output_root, segment_name);
        if exist(segment_output_dir, 'dir') ~= 7
            mkdir(segment_output_dir);
        end

        method_specs = { ...
            struct('id', 'codex_dp', 'label', 'CODEX DP', 'kind', 'curve'), ...
            struct('id', 'main_nsga2', 'label', 'NSGA-II (main)', 'kind', 'pareto'), ...
            struct('id', 'main_mopso', 'label', 'MOPSO (main)', 'kind', 'pareto'), ...
            struct('id', 'convex', 'label', 'Convex Optimisation', 'kind', 'single'), ...
            struct('id', 'milp', 'label', 'MILP', 'kind', 'single'), ...
            struct('id', 'improved_e1', 'label', 'Improved NSGA-II + improved CC_CR', 'kind', 'pareto')};

        row_structs = repmat(init_row_struct(), 0, 1);
        plot_payload = cell(0, 1);

        for method_idx = 1:numel(method_specs)
            method = method_specs{method_idx};
            switch method.id
                case 'codex_dp'
                    payload = collect_codex_payload(seg_run, target_time_s, target_tolerance_s);
                case 'main_nsga2'
                    payload = collect_main_benchmark_payload(project_root, segment_name, 'nsga2', target_time_s, target_tolerance_s);
                case 'main_mopso'
                    payload = collect_main_benchmark_payload(project_root, segment_name, 'mopso', target_time_s, target_tolerance_s);
                case 'convex'
                    payload = collect_single_result_payload(project_root, fullfile('Convex Optimisation', 'results'), segment_name, 'CO', preferred_co_solver, target_time_s);
                case 'milp'
                    payload = collect_single_result_payload(project_root, fullfile('MILP', 'results'), segment_name, 'MILP', preferred_milp_solver, target_time_s);
                case 'improved_e1'
                    payload = collect_e1_payload(project_root, segment_name, improved_config_id, target_time_s, target_tolerance_s);
                otherwise
                    payload = init_payload();
                    payload.note = sprintf('Unsupported method id: %s', method.id);
            end

            row = build_method_row(segment_name, target_time_s, method, payload);
            row_structs(end+1, 1) = row; %#ok<AGROW>

            if payload.available && ~isempty(payload.front_table)
                plot_payload{end+1, 1} = struct( ...
                    'method', method, ...
                    'payload', payload); %#ok<AGROW>
            end
        end

        row_table = struct2table(row_structs);
        row_table = add_gap_columns(row_table);

        segment_csv = fullfile(segment_output_dir, sprintf('%s_method_comparison.csv', segment_name));
        writetable(row_table, segment_csv);

        segment_png = fullfile(segment_output_dir, sprintf('%s_method_comparison.png', segment_name));
        save_segment_plot(segment_name, target_time_s, plot_payload, segment_png, showFigures);

        all_rows = [all_rows; row_table]; %#ok<AGROW>
        summary(end+1, 1) = struct( ...
            'segment', segment_name, ...
            'segment_csv', segment_csv, ...
            'segment_png', segment_png, ...
            'available_methods', sum(row_table.available), ...
            'missing_methods', strjoin(cellstr(row_table.method(~row_table.available)), ', ')); %#ok<AGROW>
    end

    aggregate_csv = fullfile(output_root, 'codex_vs_other_methods_summary.csv');
    if ~isempty(all_rows)
        writetable(all_rows, aggregate_csv);
    end
end

function payload = collect_codex_payload(seg_run, target_time_s, target_tolerance_s)
    payload = init_payload();
    payload.source_file = string(seg_run.source_file);

    target_results = seg_run.target_results;
    if isempty(target_results)
        payload.note = 'CODEX target_results is empty.';
        return;
    end

    rows = [];
    for idx = 1:numel(target_results)
        item = target_results(idx);
        if ~isfield(item, 'T_total') || ~isfield(item, 'E_total_kWh')
            continue;
        end

        actual_time = double(item.T_total);
        energy_kwh = double(item.E_total_kWh);
        if ~(isfinite(actual_time) && isfinite(energy_kwh))
            continue;
        end

        runtime_s = NaN;
        if isfield(item, 'solve_time_s')
            runtime_s = double(item.solve_time_s);
        end
        source_target_s = NaN;
        if isfield(item, 'target_T')
            source_target_s = double(item.target_T);
        end
        within_tolerance = true;
        if isfield(item, 'within_tolerance')
            within_tolerance = logical(item.within_tolerance);
        end

        rows = [rows; [actual_time, energy_kwh, runtime_s, source_target_s, double(within_tolerance)]]; %#ok<AGROW>
    end

    if isempty(rows)
        payload.note = 'CODEX target_results has no finite time-energy pairs.';
        return;
    end

    payload.front_table = array2table(rows, 'VariableNames', {'time_s', 'energy_kWh', 'runtime_s', 'source_target_s', 'is_feasible'});
    payload.front_table = sortrows(payload.front_table, {'time_s', 'energy_kWh'});
    payload.available = true;
    payload.point_rule = 'min-energy within target window; fallback closest time';
    payload.note = sprintf('CODEX total search time = %.3f s', double(seg_run.total_search_time_s));
    payload = assign_representative_point(payload, target_time_s, target_tolerance_s);
end

function payload = collect_main_benchmark_payload(project_root, segment_name, solver_name, target_time_s, target_tolerance_s)
    payload = init_payload();
    benchmark_dir = fullfile(project_root, 'main', 'results', segment_name);
    listing = dir(fullfile(benchmark_dir, 'benchmark_results*.mat'));
    if isempty(listing)
        payload.note = 'Benchmark MAT file not found.';
        return;
    end

    matfile = fullfile(benchmark_dir, listing(1).name);
    data = load(matfile, 'results');
    if ~isfield(data, 'results') || ~isstruct(data.results)
        payload.note = 'Variable results not found in benchmark MAT.';
        return;
    end

    runs = data.results(strcmp({data.results.solver}, solver_name));
    if isempty(runs)
        payload.note = sprintf('No %s runs found in benchmark MAT.', solver_name);
        return;
    end

    [front_table, point_count] = build_union_front_table_from_runs(runs);
    if isempty(front_table)
        payload.note = sprintf('No readable front points found for %s.', solver_name);
        return;
    end

    payload.available = true;
    payload.source_file = string(matfile);
    payload.front_table = front_table;
    payload.num_front_points = point_count;
    payload.point_rule = 'union Pareto front; min-energy within target window; fallback closest time';
    payload = assign_representative_point(payload, target_time_s, target_tolerance_s);
end

function payload = collect_single_result_payload(project_root, method_root_parts, segment_name, file_prefix, preferred_solver, target_time_s)
    payload = init_payload();
    result_dir = fullfile(project_root, method_root_parts, segment_name);
    if exist(result_dir, 'dir') ~= 7
        payload.note = 'Result directory not found.';
        return;
    end

    listing = dir(fullfile(result_dir, sprintf('%s_%s*_result.mat', file_prefix, segment_name)));
    if isempty(listing)
        payload.note = 'Result MAT file not found.';
        return;
    end

    selected_idx = select_preferred_file(listing, preferred_solver);
    matfile = fullfile(result_dir, listing(selected_idx).name);
    data = load(matfile, 'result');
    if ~isfield(data, 'result') || ~isstruct(data.result)
        payload.note = 'Variable result not found in MAT file.';
        return;
    end

    result = data.result;
    runtime_s = NaN;
    rows = [double(result.T_actual), double(result.E_kWh), runtime_s, double(result.T_target), 1];
    payload.available = true;
    payload.source_file = string(matfile);
    payload.front_table = array2table(rows, 'VariableNames', {'time_s', 'energy_kWh', 'runtime_s', 'source_target_s', 'is_feasible'});
    payload.point_rule = 'single deterministic solution';
    payload.representative_time_s = double(result.T_actual);
    payload.representative_energy_kWh = double(result.E_kWh);
    payload.representative_runtime_s = runtime_s;
    payload.representative_source_seed = NaN;
    payload.representative_status = string(result.status);
    payload.num_front_points = 1;
    payload.note = sprintf('Solver=%s | Target=%.3f s', string(result.solver), double(result.T_target));
    payload.target_time_s = target_time_s;
end

function payload = collect_e1_payload(project_root, segment_name, config_id, target_time_s, target_tolerance_s)
    payload = init_payload();
    matfile = fullfile(project_root, 'experiment_results', segment_name, 'E1_all_configs.mat');
    if exist(matfile, 'file') ~= 2
        payload.note = sprintf('E1_all_configs.mat not found for %s.', segment_name);
        return;
    end

    data = load(matfile, 'results');
    if ~isfield(data, 'results') || ~isstruct(data.results)
        payload.note = 'Variable results not found in E1 MAT.';
        return;
    end

    runs = data.results(strcmp({data.results.config_id}, config_id));
    if isempty(runs)
        payload.note = sprintf('Config %s not found in E1 MAT.', config_id);
        return;
    end

    [front_table, point_count] = build_union_front_table_from_runs(runs);
    if isempty(front_table)
        payload.note = sprintf('No readable front points found for E1 config %s.', config_id);
        return;
    end

    payload.available = true;
    payload.source_file = string(matfile);
    payload.front_table = front_table;
    payload.num_front_points = point_count;
    payload.point_rule = sprintf('union Pareto front from E1 config %s; min-energy within target window; fallback closest time', config_id);
    payload = assign_representative_point(payload, target_time_s, target_tolerance_s);
end

function [front_table, point_count] = build_union_front_table_from_runs(runs)
    union_rows = table();

    for idx = 1:numel(runs)
        F = extract_front(runs(idx));
        if isempty(F)
            continue;
        end

        seed = NaN;
        runtime_s = NaN;
        if isfield(runs(idx), 'seed')
            seed = double(runs(idx).seed);
        end
        if isfield(runs(idx), 'runtime')
            runtime_s = double(runs(idx).runtime);
        end

        n_rows = size(F, 1);
        chunk = table( ...
            double(F(:, 1)), ...
            double(F(:, 2)), ...
            repmat(runtime_s, n_rows, 1), ...
            repmat(seed, n_rows, 1), ...
            true(n_rows, 1), ...
            'VariableNames', {'time_s', 'energy_kWh', 'runtime_s', 'source_seed', 'is_feasible'});
        union_rows = [union_rows; chunk]; %#ok<AGROW>
    end

    if isempty(union_rows)
        front_table = table();
        point_count = 0;
        return;
    end

    F_union = [union_rows.time_s, union_rows.energy_kWh];
    keep_mask = nondominated_mask(F_union);
    front_table = union_rows(keep_mask, :);
    front_table = sortrows(front_table, {'time_s', 'energy_kWh'});
    [~, unique_idx] = unique(front_table(:, {'time_s', 'energy_kWh'}), 'rows', 'stable');
    front_table = front_table(unique_idx, :);
    point_count = height(front_table);
end

function payload = assign_representative_point(payload, target_time_s, target_tolerance_s)
    if isempty(payload.front_table)
        return;
    end

    front_table = payload.front_table;
    if ismember('is_feasible', front_table.Properties.VariableNames)
        feasible_mask = logical(front_table.is_feasible);
    else
        feasible_mask = true(height(front_table), 1);
    end
    if any(feasible_mask)
        candidate_table = front_table(feasible_mask, :);
    else
        candidate_table = front_table;
    end

    time_delta = abs(candidate_table.time_s - target_time_s);
    in_window = time_delta <= target_tolerance_s;
    if any(in_window)
        window_table = candidate_table(in_window, :);
        rel_idx = best_row_index([window_table.energy_kWh, abs(window_table.time_s - target_time_s)]);
        chosen = window_table(rel_idx, :);
    else
        rel_idx = best_row_index([time_delta, candidate_table.energy_kWh]);
        chosen = candidate_table(rel_idx, :);
    end

    payload.representative_time_s = chosen.time_s(1);
    payload.representative_energy_kWh = chosen.energy_kWh(1);
    if ismember('runtime_s', chosen.Properties.VariableNames)
        payload.representative_runtime_s = chosen.runtime_s(1);
    end
    if ismember('source_seed', chosen.Properties.VariableNames)
        payload.representative_source_seed = chosen.source_seed(1);
    end
    if ismember('source_target_s', chosen.Properties.VariableNames)
        payload.representative_source_target_s = chosen.source_target_s(1);
    end
    if abs(payload.representative_time_s - target_time_s) <= target_tolerance_s
        payload.representative_status = "within target window";
    else
        payload.representative_status = "closest to target";
    end
    payload.target_time_s = target_time_s;
    payload.num_front_points = max(payload.num_front_points, height(front_table));
end

function row = build_method_row(segment_name, target_time_s, method, payload)
    row = init_row_struct();
    row.segment = string(segment_name);
    row.method = string(method.label);
    row.available = payload.available;
    row.target_time_s = target_time_s;
    row.point_rule = string(payload.point_rule);
    row.source_file = string(payload.source_file);
    row.representative_time_s = payload.representative_time_s;
    row.representative_energy_kWh = payload.representative_energy_kWh;
    row.representative_runtime_s = payload.representative_runtime_s;
    row.source_seed = payload.representative_source_seed;
    row.num_front_points = payload.num_front_points;
    row.status = string(payload.representative_status);
    row.note = string(payload.note);

    if payload.available && ~isempty(payload.front_table)
        row.front_min_time_s = min(payload.front_table.time_s);
        row.front_max_time_s = max(payload.front_table.time_s);
        row.front_min_energy_kWh = min(payload.front_table.energy_kWh);
        row.front_max_energy_kWh = max(payload.front_table.energy_kWh);
    end
end

function row_table = add_gap_columns(row_table)
    row_table.time_gap_vs_target_s = row_table.representative_time_s - row_table.target_time_s;
    row_table.energy_gap_vs_best_pct = nan(height(row_table), 1);

    valid_energy = row_table.representative_energy_kWh(row_table.available & isfinite(row_table.representative_energy_kWh));
    if ~isempty(valid_energy)
        best_energy = min(valid_energy);
        row_table.energy_gap_vs_best_pct(row_table.available) = 100 * (row_table.representative_energy_kWh(row_table.available) - best_energy) / best_energy;
    end
end

function save_segment_plot(segment_name, target_time_s, plot_payload, output_png, showFigures)
    fig = figure('Visible', figure_visibility(showFigures), 'Color', 'w', ...
        'Name', sprintf('Method Comparison | %s', segment_name), 'NumberTitle', 'off');
    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    colors = containers.Map( ...
        {'CODEX DP', 'NSGA-II (main)', 'MOPSO (main)', 'Convex Optimisation', 'MILP', 'Improved NSGA-II + improved CC_CR'}, ...
        {[0.10 0.35 0.85], [0.00 0.55 0.30], [0.85 0.45 0.10], [0.55 0.10 0.70], [0.80 0.15 0.15], [0.15 0.60 0.70]});
    markers = containers.Map( ...
        {'CODEX DP', 'NSGA-II (main)', 'MOPSO (main)', 'Convex Optimisation', 'MILP', 'Improved NSGA-II + improved CC_CR'}, ...
        {'o', 's', 'd', '^', 'v', 'p'});

    xline(ax, target_time_s, '--k', sprintf('T_{sched}=%.2f s', target_time_s), ...
        'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.0);

    for idx = 1:numel(plot_payload)
        item = plot_payload{idx};
        label = char(item.method.label);
        payload = item.payload;
        color = colors(label);
        marker = markers(label);

        front_table = payload.front_table;
        plot(ax, front_table.time_s, front_table.energy_kWh, '-', ...
            'Color', color, 'LineWidth', 1.5, 'DisplayName', label);
        scatter(ax, front_table.time_s, front_table.energy_kWh, 28, ...
            'Marker', marker, ...
            'MarkerEdgeColor', color, ...
            'MarkerFaceColor', color, ...
            'MarkerFaceAlpha', 0.30, ...
            'DisplayName', sprintf('%s points', label));
        scatter(ax, payload.representative_time_s, payload.representative_energy_kWh, 88, ...
            'Marker', marker, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', color, ...
            'LineWidth', 1.1, ...
            'DisplayName', sprintf('%s selected', label));
    end

    xlabel(ax, 'Journey time (s)');
    ylabel(ax, 'Energy (kWh)');
    title(ax, sprintf('CODEX vs Other Methods | %s', segment_name));
    legend(ax, 'Location', 'bestoutside');

    saveas(fig, output_png);
    if ~showFigures
        close(fig);
    end
end

function bundle = load_codex_bundle(project_root)
    aggregate_file = fullfile(project_root, 'Dynamic_Programming', 'CODEX', 'results', 'all_segments_codex_results.mat');
    if exist(aggregate_file, 'file') ~= 2
        error('CODEX aggregate file not found: %s', aggregate_file);
    end

    data = load(aggregate_file, 'segment_runs');
    if ~isfield(data, 'segment_runs') || ~isstruct(data.segment_runs)
        error('segment_runs not found in %s', aggregate_file);
    end

    segment_runs = data.segment_runs;
    for idx = 1:numel(segment_runs)
        segment_runs(idx).source_file = aggregate_file;
        if ~isfield(segment_runs(idx), 'total_search_time_s')
            segment_runs(idx).total_search_time_s = NaN;
        end
    end

    bundle = struct();
    bundle.aggregate_file = aggregate_file;
    bundle.segment_runs = segment_runs;
end

function selected_runs = resolve_selected_runs(segment_runs, selected_segments)
    if ischar(selected_segments) || (isstring(selected_segments) && isscalar(selected_segments))
        token = char(string(selected_segments));
        if strcmpi(token, 'auto') || strcmpi(token, 'all')
            selected_runs = segment_runs;
            return;
        end
        selected_segments = {token};
    elseif isstring(selected_segments)
        selected_segments = cellstr(selected_segments(:));
    elseif ~iscell(selected_segments)
        error('selected_segments must be ''auto'', ''all'', a string, or a cell array of strings.');
    end

    available_names = {segment_runs.segment};
    available_names = cellfun(@(item) char(string(item.name)), available_names, 'UniformOutput', false);
    keep_mask = false(1, numel(segment_runs));
    for idx = 1:numel(selected_segments)
        match_idx = find(strcmpi(available_names, char(string(selected_segments{idx}))), 1);
        if isempty(match_idx)
            error('Requested segment %s is not available in CODEX aggregate results.', char(string(selected_segments{idx})));
        end
        keep_mask(match_idx) = true;
    end
    selected_runs = segment_runs(keep_mask);
end

function F = extract_front(run)
    F = [];
    if isfield(run, 'F') && ~isempty(run.F)
        F = double(run.F);
        return;
    end
    if isfield(run, 'AllF') && ~isempty(run.AllF)
        F = double(run.AllF);
    end
end

function keep_mask = nondominated_mask(F)
    n = size(F, 1);
    keep_mask = true(n, 1);
    if isempty(F)
        return;
    end

    invalid_mask = any(~isfinite(F), 2);
    keep_mask(invalid_mask) = false;
    for i = 1:n
        if ~keep_mask(i)
            continue;
        end
        dominates_i = all(F <= F(i, :), 2) & any(F < F(i, :), 2);
        dominates_i(i) = false;
        if any(dominates_i)
            keep_mask(i) = false;
        end
    end
end

function selected_idx = select_preferred_file(listing, preferred_solver)
    selected_idx = 1;
    if isempty(preferred_solver)
        return;
    end

    solver_token = lower(string(preferred_solver));
    names = string({listing.name});
    match_idx = find(contains(lower(names), "_" + solver_token + "_"), 1);
    if isempty(match_idx)
        match_idx = find(contains(lower(names), solver_token), 1);
    end
    if ~isempty(match_idx)
        selected_idx = match_idx;
    end
end

function visibility = figure_visibility(showFigures)
    if showFigures
        visibility = 'on';
    else
        visibility = 'off';
    end
end

function idx = best_row_index(values)
    [~, order] = sortrows(values);
    idx = order(1);
end

function payload = init_payload()
    payload = struct( ...
        'available', false, ...
        'source_file', "", ...
        'front_table', table(), ...
        'representative_time_s', NaN, ...
        'representative_energy_kWh', NaN, ...
        'representative_runtime_s', NaN, ...
        'representative_source_seed', NaN, ...
        'representative_source_target_s', NaN, ...
        'representative_status', "missing", ...
        'point_rule', "", ...
        'num_front_points', 0, ...
        'note', "", ...
        'target_time_s', NaN);
end

function row = init_row_struct()
    row = struct( ...
        'segment', "", ...
        'method', "", ...
        'available', false, ...
        'target_time_s', NaN, ...
        'representative_time_s', NaN, ...
        'representative_energy_kWh', NaN, ...
        'representative_runtime_s', NaN, ...
        'time_gap_vs_target_s', NaN, ...
        'energy_gap_vs_best_pct', NaN, ...
        'source_seed', NaN, ...
        'num_front_points', 0, ...
        'front_min_time_s', NaN, ...
        'front_max_time_s', NaN, ...
        'front_min_energy_kWh', NaN, ...
        'front_max_energy_kWh', NaN, ...
        'status', "missing", ...
        'point_rule', "", ...
        'source_file', "", ...
        'note', "");
end

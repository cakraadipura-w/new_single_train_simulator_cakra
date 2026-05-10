%% read_benchmark_results.m
% Run this file directly from MATLAB Editor.
%
% Quick settings:
%   selected_results = {'IS01'};
%   selected_results = {'IS01','IS02'};
%   selected_results = 'all';
%   selected_results = 'IS01';
%   selected_seeds = 'all';
%   selected_seeds = [1 2];
%
% Notes:
%   - Per-route CSV/PNG stay inside each route folder under results/.
%   - Aggregate CSV is written to outputRoot.
clear

if ~exist('selected_results', 'var') || isempty(selected_results)
    selected_results = {'IS06'};
end

this_dir = fileparts(mfilename('fullpath'));
default_root = resolve_default_root(this_dir);

if ~exist('outputRoot', 'var') || isempty(outputRoot)
    outputRoot = default_root;
end

if ~exist('showFigures', 'var') || isempty(showFigures)
    showFigures = true;
end

if ~exist('selected_seeds', 'var') || isempty(selected_seeds)
    selected_seeds = 'all';
end

if ~exist('enableParetoControls', 'var') || isempty(enableParetoControls)
    enableParetoControls = true;
end

summary = read_benchmark_results_impl(selected_results, outputRoot, showFigures, selected_seeds, enableParetoControls);

if ~isempty(summary)
    summary_table = struct2table(summary);
    disp(summary_table(:, {'route', 'num_runs', 'summary_csv'}));
end

function summary = read_benchmark_results_impl(matInput, outputRoot, showFigures, selectedSeeds, enableParetoControls)
%READ_BENCHMARK_RESULTS_IMPL Summarize benchmark MAT files into CSV and plots.

    this_dir = fileparts(mfilename('fullpath'));
    default_root = resolve_default_root(this_dir);

    if nargin < 1 || isempty(matInput)
        matInput = default_root;
    end
    if nargin < 2 || isempty(outputRoot)
        outputRoot = default_root;
    end
    if nargin < 4 || isempty(selectedSeeds)
        selectedSeeds = 'all';
    end
    if nargin < 5 || isempty(enableParetoControls)
        enableParetoControls = logical(showFigures);
    end

    if ~exist(outputRoot, 'dir')
        mkdir(outputRoot);
    end

    files = resolve_mat_files(matInput, default_root);
    if isempty(files)
        error('No benchmark MAT files found for input: %s', stringify_input(matInput));
    end
    if nargin < 3 || isempty(showFigures)
        showFigures = numel(files) <= 2;
    end

    summary = repmat(struct( ...
        'matfile', '', ...
        'route', '', ...
        'summary_csv', '', ...
        'pareto_png', '', ...
        'metrics_png', '', ...
        'num_runs', 0), 0, 1);

    all_rows = table();

    for file_idx = 1:numel(files)
        matfile = files{file_idx};
        payload = load(matfile, 'results');
        if ~isfield(payload, 'results') || ~isstruct(payload.results)
            warning('Skipping %s because variable "results" was not found.', matfile);
            continue;
        end

        source_run_indices = (1:numel(payload.results))';
        runs = payload.results;
        if isempty(runs)
            warning('Skipping %s because results is empty.', matfile);
            continue;
        end

        [runs, source_run_indices] = filter_runs_by_seed(runs, source_run_indices, selectedSeeds);
        if isempty(runs)
            warning('Skipping %s because no runs matched selected_seeds=%s.', matfile, stringify_input(selectedSeeds));
            continue;
        end

        [route_name, base_name] = infer_labels(matfile);
        row_table = build_summary_table(runs, route_name, matfile, source_run_indices);
        if isempty(row_table)
            warning('Skipping %s because no readable runs were found.', matfile);
            continue;
        end

        summary_csv = fullfile(fileparts(matfile), [base_name '_summary.csv']);
        pareto_png = fullfile(fileparts(matfile), [base_name '_pareto.png']);
        metrics_png = fullfile(fileparts(matfile), [base_name '_metrics.png']);

        writetable(row_table, summary_csv);
        save_pareto_plot(runs, source_run_indices, route_name, base_name, pareto_png, showFigures, enableParetoControls);
        save_metric_plot(row_table, route_name, base_name, metrics_png, showFigures);

        if isempty(all_rows)
            all_rows = row_table;
        else
            all_rows = [all_rows; row_table]; %#ok<AGROW>
        end

        summary(end+1, 1) = struct( ...
            'matfile', matfile, ...
            'route', route_name, ...
            'summary_csv', summary_csv, ...
            'pareto_png', pareto_png, ...
            'metrics_png', metrics_png, ...
            'num_runs', height(row_table)); %#ok<AGROW>

        fprintf('Processed %s -> %s\n', matfile, summary_csv);
    end

    if isempty(summary)
        error('No readable benchmark MAT files were processed.');
    end

    aggregate_csv = fullfile(outputRoot, 'benchmark_results_summary_all.csv');
    writetable(all_rows, aggregate_csv);
    fprintf('Saved aggregate summary: %s\n', aggregate_csv);
end

function default_root = resolve_default_root(this_dir)
    [~, folder_name] = fileparts(this_dir);
    if strcmpi(folder_name, 'results')
        default_root = this_dir;
    else
        default_root = fullfile(this_dir, 'results');
    end
end

function files = resolve_mat_files(matInput, default_root)
    files = {};

    if iscell(matInput)
        for idx = 1:numel(matInput)
            files = [files; resolve_mat_files(matInput{idx}, default_root)]; %#ok<AGROW>
        end
        files = unique(files, 'stable');
        return;
    end

    input_text = char(string(matInput));
    if strcmpi(input_text, 'all')
        files = collect_benchmark_files(default_root);
        return;
    end

    if isfolder(input_text)
        files = collect_benchmark_files(input_text);
        return;
    end

    candidate_dir = fullfile(default_root, input_text);
    if isfolder(candidate_dir)
        files = collect_benchmark_files(candidate_dir);
        return;
    end

    if isfile(input_text)
        files = {input_text};
        return;
    end

    candidate_file = fullfile(default_root, input_text);
    if isfile(candidate_file)
        files = {candidate_file};
        return;
    end

    match = dir(input_text);
    if isempty(match)
        match = dir(fullfile(default_root, input_text));
    end
    for idx = 1:numel(match)
        if ~match(idx).isdir
            files{end+1, 1} = fullfile(match(idx).folder, match(idx).name); %#ok<AGROW>
        end
    end
    files = unique(files, 'stable');
end

function files = collect_benchmark_files(root_dir)
    entries = dir(root_dir);
    files = {};

    for idx = 1:numel(entries)
        name = entries(idx).name;
        if entries(idx).isdir
            if strcmp(name, '.') || strcmp(name, '..')
                continue;
            end
            sub_files = collect_benchmark_files(fullfile(root_dir, name));
            if ~isempty(sub_files)
                files = [files; sub_files]; %#ok<AGROW>
            end
        elseif startsWith(name, 'benchmark_results_') && endsWith(name, '.mat')
            files{end+1, 1} = fullfile(root_dir, name); %#ok<AGROW>
        end
    end
end

function row_table = build_summary_table(runs, route_name, matfile, source_run_indices)
    n = numel(runs);

    route_col = repmat(string(route_name), n, 1);
    matfile_col = repmat(string(matfile), n, 1);
    solver_col = strings(n, 1);
    strategy_col = strings(n, 1);
    seed_col = nan(n, 1);
    pop_col = nan(n, 1);
    dim_col = nan(n, 1);
    iter_col = nan(n, 1);
    runtime_col = nan(n, 1);
    hv_col = nan(n, 1);
    igd_col = nan(n, 1);
    nf1_col = nan(n, 1);
    best_time_col = nan(n, 1);
    best_energy_col = nan(n, 1);
    t_at_best_energy_col = nan(n, 1);

    for idx = 1:n
        run = runs(idx);

        if isfield(run, 'solver')
            solver_col(idx) = string(run.solver);
        else
            solver_col(idx) = "unknown";
        end
        if isfield(run, 'strategy')
            strategy_col(idx) = string(run.strategy);
        else
            strategy_col(idx) = "unknown";
        end
        if isfield(run, 'seed')
            seed_col(idx) = double(run.seed);
        end
        if isfield(run, 'pop_size')
            pop_col(idx) = double(run.pop_size);
        end
        if isfield(run, 'dim')
            dim_col(idx) = double(run.dim);
        end
        if isfield(run, 'iterations')
            iter_col(idx) = double(run.iterations);
        end
        if isfield(run, 'runtime')
            runtime_col(idx) = double(run.runtime);
        end
        if isfield(run, 'HV')
            hv_col(idx) = double(run.HV);
        end
        if isfield(run, 'IGD')
            igd_col(idx) = double(run.IGD);
        end
        if isfield(run, 'Nf1')
            nf1_col(idx) = double(run.Nf1);
        end

        F = extract_front(run);
        if ~isempty(F)
            best_time_col(idx) = min(F(:, 1));
            [best_energy_col(idx), best_idx] = min(F(:, 2));
            t_at_best_energy_col(idx) = F(best_idx, 1);
        end
    end

    row_table = table( ...
        route_col, ...
        matfile_col, ...
        source_run_indices(:), ...
        solver_col, ...
        strategy_col, ...
        seed_col, ...
        pop_col, ...
        dim_col, ...
        iter_col, ...
        runtime_col, ...
        hv_col, ...
        igd_col, ...
        nf1_col, ...
        best_time_col, ...
        best_energy_col, ...
        t_at_best_energy_col, ...
        'VariableNames', { ...
            'route', 'matfile', 'run_index', 'solver', 'strategy', 'seed', ...
            'pop_size', 'dim', 'iterations', 'runtime_s', 'HV', 'IGD', 'Nf1', ...
            'best_time_s', 'best_energy_kWh', 'time_at_best_energy_s'});
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

function save_pareto_plot(runs, source_run_indices, route_name, base_name, output_png, showFigures, enableParetoControls)
    fig = figure('Visible', figure_visibility(showFigures), 'Color', 'w', ...
        'Name', sprintf('Pareto | %s', route_name), 'NumberTitle', 'off');
    if showFigures && enableParetoControls
        ax = axes(fig, 'Position', [0.08 0.12 0.60 0.80]);
    else
        ax = axes(fig);
    end
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    solvers = unique(string({runs.solver}), 'stable');
    colors = lines(max(numel(solvers), 1));
    seed_values = extract_seed_values(runs);
    markers = {'o','s','d','^','v','>','<','p','h','x','+','*'};
    plotted = false;
    plot_handles = gobjects(0);
    plot_meta = struct('solver', {}, 'seed', {}, 'runIndex', {}, 'seriesLabel', {});

    for idx = 1:numel(runs)
        run = runs(idx);
        F = extract_front(run);
        if isempty(F)
            continue;
        end

        solver = extract_solver_name(run);
        seed = extract_seed_number(run);
        solver_idx = find(solvers == solver, 1);
        if isempty(solver_idx)
            solver_idx = 1;
        end
        seed_idx = find(seed_values == seed, 1);
        if isempty(seed_idx)
            seed_idx = 1;
        end

        series_label = build_series_label(solver, seed, source_run_indices(idx));
        h = scatter(ax, F(:, 1), F(:, 2), 28, ...
            'Marker', markers{mod(seed_idx - 1, numel(markers)) + 1}, ...
            'MarkerFaceColor', colors(solver_idx, :), ...
            'MarkerEdgeColor', colors(solver_idx, :), ...
            'MarkerFaceAlpha', 0.45, ...
            'DisplayName', series_label);
        h.UserData = struct( ...
            'solver', char(solver), ...
            'seed', seed, ...
            'runIndex', source_run_indices(idx), ...
            'seriesLabel', series_label);

        plot_handles(end+1, 1) = h; %#ok<AGROW>
        plot_meta(end+1, 1) = h.UserData; %#ok<AGROW>
        plotted = true;
    end

    xlabel(ax, 'Running time (s)');
    ylabel(ax, 'Energy (kWh)');
    title(ax, sprintf('%s | %s', route_name, strrep(base_name, '_', '\_')));
    if plotted
        lgd = legend(ax, plot_handles, 'Location', 'eastoutside', 'Interpreter', 'none');
        saveas(fig, output_png);
        if showFigures
            try
                dcm = datacursormode(fig);
                set(dcm, 'Enable', 'on', 'UpdateFcn', @pareto_datatip);
            catch
                % Keep figure usable even if custom data tips are unsupported.
            end
        end

        if showFigures && enableParetoControls
            if isgraphics(lgd)
                delete(lgd);
            end
            add_toggle_controls();
        end

        if showFigures
            figure(fig);
            drawnow;
        end
    else
        close(fig);
        return;
    end
    if ~showFigures
        close(fig);
    end

    function add_toggle_controls()
        panel = uipanel('Parent', fig, 'Title', 'On / Off Compare', ...
            'Units', 'normalized', 'Position', [0.70 0.06 0.28 0.88]);

        uicontrol(panel, 'Style', 'pushbutton', 'String', 'All On', ...
            'Units', 'normalized', 'Position', [0.05 0.94 0.40 0.05], ...
            'Callback', @(~,~) set_all_series(true));
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'All Off', ...
            'Units', 'normalized', 'Position', [0.55 0.94 0.40 0.05], ...
            'Callback', @(~,~) set_all_series(false));

        uicontrol(panel, 'Style', 'text', 'String', 'Show only seed:', ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
            'Position', [0.05 0.88 0.90 0.04]);

        unique_seeds = unique([plot_meta.seed], 'stable');
        y_seed = 0.83;
        for seed_button_idx = 1:numel(unique_seeds)
            seed_value = unique_seeds(seed_button_idx);
            uicontrol(panel, 'Style', 'pushbutton', ...
                'String', sprintf('Seed %d', seed_value), ...
                'Units', 'normalized', 'Position', [0.05 y_seed 0.27 0.045], ...
                'Callback', @(~,~) show_only_seed(seed_value));
            y_seed = y_seed - 0.055;
        end

        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Show All Seeds', ...
            'Units', 'normalized', 'Position', [0.36 0.83 0.59 0.045], ...
            'Callback', @(~,~) set_all_series(true));

        uicontrol(panel, 'Style', 'text', 'String', 'Runs:', ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
            'Position', [0.05 0.62 0.90 0.04]);

        n_series = numel(plot_handles);
        series_height = min(0.048, 0.56 / max(n_series, 1));
        y_series = 0.58;
        series_boxes = gobjects(n_series, 1);

        for box_idx = 1:n_series
            series_boxes(box_idx) = uicontrol(panel, 'Style', 'checkbox', ...
                'String', plot_meta(box_idx).seriesLabel, ...
                'Units', 'normalized', 'Value', 1, ...
                'HorizontalAlignment', 'left', ...
                'Position', [0.05 y_series 0.90 series_height], ...
                'Callback', @(src,~) set(plot_handles(box_idx), 'Visible', on_off(src.Value)));
            y_series = y_series - series_height - 0.008;
            if y_series < 0.02
                break;
            end
        end

        function set_all_series(is_visible)
            for control_idx = 1:numel(series_boxes)
                if isgraphics(series_boxes(control_idx))
                    series_boxes(control_idx).Value = double(is_visible);
                end
                if isgraphics(plot_handles(control_idx))
                    set(plot_handles(control_idx), 'Visible', on_off(is_visible));
                end
            end
        end

        function show_only_seed(seed_value)
            for control_idx = 1:numel(series_boxes)
                is_visible = plot_meta(control_idx).seed == seed_value;
                if isgraphics(series_boxes(control_idx))
                    series_boxes(control_idx).Value = double(is_visible);
                end
                if isgraphics(plot_handles(control_idx))
                    set(plot_handles(control_idx), 'Visible', on_off(is_visible));
                end
            end
        end
    end

    function txt = pareto_datatip(~, event)
        meta = event.Target.UserData;
        pos = event.Position;
        txt = { ...
            sprintf('Solver: %s', meta.solver), ...
            sprintf('Seed: %d', meta.seed), ...
            sprintf('Run: %d', meta.runIndex), ...
            sprintf('Time: %.6f s', pos(1)), ...
            sprintf('Energy: %.6f kWh', pos(2))};
    end
end

function save_metric_plot(rows, route_name, base_name, output_png, showFigures)
    fig = figure('Visible', figure_visibility(showFigures), 'Color', 'w', ...
        'Name', sprintf('Metrics | %s', route_name), 'NumberTitle', 'off');

    metric_names = {'HV', 'IGD', 'runtime_s', 'Nf1'};
    metric_titles = {'HV', 'IGD', 'Runtime (s)', 'Pareto count'};

    for idx = 1:numel(metric_names)
        subplot(2, 2, idx);
        [solvers, means, stds] = summarize_metric(rows, metric_names{idx});
        if isempty(solvers)
            axis off;
            title(sprintf('%s (no data)', metric_titles{idx}));
            continue;
        end

        x = 1:numel(solvers);
        bar(x, means, 'FaceColor', [0.25 0.45 0.75]); hold on;
        errorbar(x, means, stds, '.', 'Color', [0 0 0], 'LineWidth', 1);
        set(gca, 'XTick', x, 'XTickLabel', cellstr(solvers));
        xtickangle(20);
        title(metric_titles{idx});
        grid on;
    end

    seed_values = unique(rows.seed(isfinite(rows.seed)));
    sgtitle(sprintf('%s | %s | seeds=%s', route_name, strrep(base_name, '_', '\_'), seed_values_text(seed_values)));
    saveas(fig, output_png);
    if showFigures
        figure(fig);
        drawnow;
    else
        close(fig);
    end
end

function [solvers, means, stds] = summarize_metric(rows, metric_name)
    solver_values = string(rows.solver);
    solvers = unique(solver_values, 'stable');
    means = nan(numel(solvers), 1);
    stds = nan(numel(solvers), 1);

    values = rows.(metric_name);
    keep_mask = false(numel(solvers), 1);

    for idx = 1:numel(solvers)
        mask = solver_values == solvers(idx) & isfinite(values);
        if ~any(mask)
            continue;
        end
        data = values(mask);
        means(idx) = mean(data);
        stds(idx) = std(data);
        keep_mask(idx) = true;
    end

    solvers = solvers(keep_mask);
    means = means(keep_mask);
    stds = stds(keep_mask);
end

function [route_name, base_name] = infer_labels(matfile)
    [folder_path, base_name, ~] = fileparts(matfile);
    [~, route_name] = fileparts(folder_path);
    if isempty(route_name)
        route_name = 'results';
    end
end

function text_value = stringify_input(input_value)
    if iscell(input_value)
        text_value = strjoin(cellfun(@char, cellstr(string(input_value)), 'UniformOutput', false), ', ');
    else
        text_value = char(string(input_value));
    end
end

function [runs_filtered, source_run_indices] = filter_runs_by_seed(runs, source_run_indices, selectedSeeds)
    runs_filtered = runs;
    if nargin < 3 || isempty(selectedSeeds)
        return;
    end

    if (ischar(selectedSeeds) || (isstring(selectedSeeds) && isscalar(selectedSeeds))) && strcmpi(string(selectedSeeds), 'all')
        return;
    end

    seed_values = normalize_seed_selector(selectedSeeds);
    keep_mask = false(size(runs));
    for idx = 1:numel(runs)
        if isfield(runs(idx), 'seed') && ~isempty(runs(idx).seed)
            keep_mask(idx) = any(double(runs(idx).seed) == seed_values);
        end
    end

    runs_filtered = runs(keep_mask);
    source_run_indices = source_run_indices(keep_mask);
end

function seed_values = normalize_seed_selector(selectedSeeds)
    if isnumeric(selectedSeeds)
        seed_values = selectedSeeds(:)';
    elseif iscell(selectedSeeds)
        seed_values = zeros(1, numel(selectedSeeds));
        for idx = 1:numel(selectedSeeds)
            item = selectedSeeds{idx};
            if isnumeric(item)
                seed_values(idx) = double(item);
            else
                seed_values(idx) = str2double(string(item));
            end
        end
    elseif isstring(selectedSeeds)
        seed_values = str2double(selectedSeeds(:))';
    elseif ischar(selectedSeeds)
        seed_values = str2double(string(selectedSeeds));
    else
        error('Unsupported selected_seeds type: %s', class(selectedSeeds));
    end

    if any(~isfinite(seed_values)) || any(seed_values ~= floor(seed_values))
        error('selected_seeds must be ''all'' or integer seed values, e.g. [1 2 3].');
    end
    seed_values = unique(seed_values, 'stable');
end

function solver = extract_solver_name(run)
    if isfield(run, 'solver') && ~isempty(run.solver)
        solver = string(run.solver);
    else
        solver = "unknown";
    end
end

function seed = extract_seed_number(run)
    if isfield(run, 'seed') && ~isempty(run.seed)
        seed = double(run.seed);
    else
        seed = NaN;
    end
end

function seed_values = extract_seed_values(runs)
    seed_values = zeros(0, 1);
    for idx = 1:numel(runs)
        seed = extract_seed_number(runs(idx));
        if isfinite(seed)
            seed_values(end+1, 1) = seed; %#ok<AGROW>
        end
    end
    seed_values = unique(seed_values, 'stable');
    if isempty(seed_values)
        seed_values = NaN;
    end
end

function label = build_series_label(solver, seed, run_index)
    if isfinite(seed)
        label = sprintf('%s | seed %d | run %d', char(solver), seed, run_index);
    else
        label = sprintf('%s | run %d', char(solver), run_index);
    end
end

function text_value = seed_values_text(seed_values)
    if isempty(seed_values)
        text_value = 'all';
        return;
    end
    text_value = strjoin(arrayfun(@(x) sprintf('%d', x), seed_values, 'UniformOutput', false), ',');
end

function state = on_off(is_visible)
    if is_visible
        state = 'on';
    else
        state = 'off';
    end
end

function visibility = figure_visibility(showFigures)
    if showFigures
        visibility = 'on';
    else
        visibility = 'off';
    end
end
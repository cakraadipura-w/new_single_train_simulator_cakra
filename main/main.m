%% main.m  (thin orchestrator)
clc; clear; clear global;
main_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(main_dir);
addpath(genpath(project_root));

% ========================= USER SETTINGS =========================
cfg = struct();

% Save all outputs under main/results, regardless of the active working folder.
cfg.output_root_dir = fullfile(main_dir, 'results');
if ~exist(cfg.output_root_dir, 'dir')
    mkdir(cfg.output_root_dir);
end

% Core data / model choices
cfg.route_direction = 'up';   % 'up' | 'down' for Guangzhou Line 7 segments
route_catalog = build_route_catalog(cfg.route_direction);

% Route selection:
%   'all'                     -> run all routes in route_catalog
%   'IS01'                    -> run a single route by key
%   {'IS01','IS03','IS05'}    -> run a subset by keys
%   [3 4 5]                   -> run a subset by route_catalog index
%   {'line7_full','IS02'}     -> mix keys / filenames from route_catalog
cfg.route_selector = {'IS06','IS07','IS08'};

cfg.rollingstock_file = 'rollingstock_Guangzhou_L7.m'; 
%cfg.rollingstock_file = 'rollingstock_Guangzhou_L7_res_upd.m'; 
cfg.driving_strategy  = "CC_CR";                      % "CC" | "CR" | "CC_CR"

% Parallel
cfg.parallel_use   = true;     % true/false
cfg.leave_one_core = false;    % true -> use (NumWorkers-1) to keep UI responsive

% Strategy switch (single-run mode)
% - false => base (simulation_fun_CC_CR_base + decision_var_NO_base)
% - true  => improved (simulation_fun_CC_CR_improved + decision_var_NO_improved)
cfg.use_improved = false;

% Optim settings (default)
cfg.pop_size   = 200;
cfg.iterations = 300;

% Run mode
cfg.run_benchmark = true;      % true -> run benchmark, false -> run one optimizer
cfg.optimizer     = "nsga2";   % "nsga2" | "mopso" | "spea2" | "moead" | "dp" | "de_moea" | "dqn" (used when run_benchmark=false)

% NSGA-II variant selection (used when optimizer is "nsga2")
% Options: 'original' | 'rl_sde' | 'improved_regular' | 'bayes_rl'
cfg.nsga2_variant = 'original';  % Default: RL+SDE hybrid (Wen+Zhang 2025)

% Benchmark settings (used when run_benchmark=true)
% NOTE: this will run EACH solver across EACH seed across EACH strategy.
cfg.bench_solvers    = {'nsga2', 'mopso'};           % e.g. {'nsga2','mopso','spea2','moead','dp','de_moea','dqn'}
cfg.bench_seeds      = 1:3;                % repeat count (stochastic)
cfg.bench_strategies = {'base'}; % compare both strategies in one benchmark
cfg.make_plots       = true;               % plots: Pareto overlay + boxplots

global time_obj_max
time_obj_max = 250;   % maximum running time (objective) limit: 250 seconds

% Fair budget helper:
% If you keep cfg.iterations for NSGA-II, the default benchmark uses ~half iterations for others.
% You can override like this:
% cfg.bench_iters = struct('nsga2',300,'mopso',150,'spea2',150,'moead',150);

% Progress log inside calculate_pop (optional)
global show_progress
show_progress = false;

% Optional sim options (forward compatible)
cfg.sim = struct();

% Global mode flag (used by simulator + decision_var dispatcher)
global use_improved
use_improved = logical(cfg.use_improved);

selected_routes = resolve_route_selection(cfg.route_selector, route_catalog);
fprintf('Selected %d route(s): %s\n', numel(selected_routes), strjoin({selected_routes.key}, ', '));

for route_idx = 1:numel(selected_routes)
    cfg.route_file = selected_routes(route_idx).file;
    route_key = selected_routes(route_idx).key;
    cfg.output_dir = fullfile(cfg.output_root_dir, route_key);
    if ~exist(cfg.output_dir, 'dir')
        mkdir(cfg.output_dir);
    end

    fprintf('\n============================================================\n');
    fprintf('Route %d/%d | %s | %s\n', route_idx, numel(selected_routes), route_key, cfg.route_file);
    fprintf('Output      | %s\n', cfg.output_dir);
    fprintf('============================================================\n');

    % ========================= SETUP =========================
    info = setup_project(cfg);

    % ========================= PARALLEL SETUP =========================
    if cfg.parallel_use
        setup_parallel_pool(cfg, info);
    end

    % ========================= RUN =========================
    % (pass resolved paths into cfg so benchmark can refresh workers when switching strategy)
    cfg.routePath        = info.routePath;
    cfg.rollingstockPath = info.rollingstockPath;
    cfg.project_root     = info.project_root;

    if cfg.run_benchmark
        results = run_benchmark_compare(cfg);

        % Generate smart filename: benchmark_results_<solvers>_<variant>_<strategies>.mat
        solvers = string(cfg.bench_solvers);
        strategies = string(cfg.bench_strategies);
        solver_str = strjoin(cellstr(solvers), '-');
        strategy_str = strjoin(strategies, '-');

        % Add variant info only when benchmarking a single NSGA-II solver
        if isscalar(solvers) && strcmp(solvers(1), 'nsga2')
            variant_str = char(cfg.nsga2_variant);
            benchmark_fname = sprintf('benchmark_results_%s_%s_%s.mat', ...
                solver_str, variant_str, strategy_str);
        else
            benchmark_fname = sprintf('benchmark_results_%s_%s.mat', ...
                solver_str, strategy_str);
        end

        benchmark_path = fullfile(cfg.output_dir, benchmark_fname);
        save(benchmark_path, 'results');
        fprintf('Saved: %s\n', benchmark_path);
    else
        OUT = run_single(cfg);

        % Generate smart filename based on optimizer + variant + strategy
        segment = extract_segment_name(cfg.route_file);
        result_fname = sprintf('result_%s_%s_%s%s.mat', ...
            OUT.optimizer, OUT.nsga2_variant, OUT.strategy, segment);

        result_path = fullfile(cfg.output_dir, result_fname);
        save(result_path, 'OUT');
        fprintf('Saved: %s\n', result_path);
    end
end

% Helper: Extract segment name from route file
function segment = extract_segment_name(route_file)
    % Extract segment identifier (e.g., '_IS01' from 'Guangzhou_Line7_IS01_0.000-1.120km.mat')
    if contains(route_file, 'IS')
        match = regexp(route_file, '_IS\d+', 'match');
        if ~isempty(match)
            segment = match{1};  % e.g., '_IS01'
        else
            segment = '';
        end
    else
        segment = '';
    end
end

function selected_routes = resolve_route_selection(route_selector, route_catalog)
    if isnumeric(route_selector)
        route_idx = unique(route_selector(:)', 'stable');
        if isempty(route_idx) || any(route_idx < 1) || any(route_idx > numel(route_catalog))
            error('route_selector index out of range. Valid index: 1..%d', numel(route_catalog));
        end
        selected_routes = route_catalog(route_idx);
        return;
    end

    if ischar(route_selector) || (isstring(route_selector) && isscalar(route_selector))
        selector_tokens = cellstr(string(route_selector));
    elseif isstring(route_selector)
        selector_tokens = cellstr(route_selector(:));
    elseif iscell(route_selector)
        selector_tokens = route_selector;
    else
        error('Unsupported route_selector type: %s', class(route_selector));
    end

    if isscalar(selector_tokens) && strcmpi(string(selector_tokens{1}), 'all')
        selected_routes = route_catalog;
        return;
    end

    route_idx = zeros(1, numel(selector_tokens));
    catalog_keys = {route_catalog.key};
    catalog_files = {route_catalog.file};

    for i = 1:numel(selector_tokens)
        token = string(selector_tokens{i});
        if strlength(token) == 0
            error('route_selector contains an empty entry.');
        end

        idx = find(strcmpi(catalog_keys, char(token)) | strcmpi(catalog_files, char(token)), 1);
        if isempty(idx)
            error('Unknown route selector: %s', char(token));
        end

        route_idx(i) = idx;
    end

    selected_routes = route_catalog(unique(route_idx, 'stable'));
end

function route_catalog = build_route_catalog(route_direction)
    line7_segments = get_guangzhou_line7_catalog(route_direction);
    route_catalog = struct( ...
        'key',  [{'xiau_liu', 'line7_full'}, {line7_segments.name}], ...
        'file', [{'xiau_liu_route.mat', 'Guangzhou_Line7.mat'}, {line7_segments.file}]);
end

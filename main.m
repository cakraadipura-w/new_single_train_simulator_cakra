%% main.m  (thin orchestrator)
clc; clear; clear global;
addpath(genpath(pwd));

% ========================= USER SETTINGS =========================
cfg = struct();

% Core data / model choices
%cfg.route_file        = 'xiau_liu_route.mat';
%cfg.rollingstock_file = 'rollingstock_Xiau_liu2.m';  % change here to switch rollingstock

%cfg.route_file        = 'Guangzhou_Line7.mat';
cfg.route_file        = 'Guangzhou_Line7_IS01_0.000-1.120km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS02_1.120-3.028km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS03_3.028-5.200km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS04_5.200-6.842km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS05_6.842-8.958km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS06_8.958-11.323km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS07_11.323-13.729km.mat';
%cfg.route_file        = 'Guangzhou_Line7_IS08_13.729-17.507km.mat';
cfg.rollingstock_file = 'rollingstock_Guangzhou_L7.m'; 
%cfg.rollingstock_file = 'rollingstock_Guangzhou_L7_res_upd.m'; 
cfg.driving_strategy  = "CC_CR";                      % "CC" | "CR" | "CC_CR"

% Parallel
cfg.parallel_use   = true;     % true/false
cfg.leave_one_core = false;    % true -> use (NumWorkers-1) to keep UI responsive

% Strategy switch (single-run mode)
% - false => base (simulation_fun_CC_CR_base + decision_var_NO_base)
% - true  => improved (simulation_fun_CC_CR_improved + decision_var_NO_improved)
cfg.use_improved = true;

% Optim settings (default)
cfg.pop_size   = 200;
cfg.iterations = 300;

% Run mode
cfg.run_benchmark = true;      % true -> run benchmark, false -> run one optimizer
cfg.optimizer     = "nsga2";   % "nsga2" | "mopso" | "spea2" | "moead" (used when run_benchmark=false)

% NSGA-II variant selection (used when optimizer is "nsga2")
% Options: 'original' | 'rl_sde' | 'improved_regular' | 'bayes_rl'
cfg.nsga2_variant = 'original';  % Default: RL+SDE hybrid (Wen+Zhang 2025)

% Benchmark settings (used when run_benchmark=true)
% NOTE: this will run EACH solver across EACH seed across EACH strategy.
cfg.bench_solvers    = {'nsga2'};           % e.g. {'nsga2','mopso','spea2','moead'}
cfg.bench_seeds      = 1:1;                % repeat count (stochastic)
cfg.bench_strategies = {'base','improved'}; % compare both strategies in one benchmark
cfg.make_plots       = true;               % plots: Pareto overlay + boxplots

global time_obj_max
time_obj_max = 250;   % maximum running time (objective) limit: 250 seconds

% Fair budget helper:
% If you keep cfg.iterations for NSGA-II, the default benchmark uses ~half iterations for others.
% You can override like this:
% cfg.bench_iters = struct('nsga2',300,'mopso',150,'spea2',150,'moead',150);

% Progress log inside calculate_pop (optional)
global show_progress
show_progress = true;

% Optional sim options (forward compatible)
cfg.sim = struct();

% Global mode flag (used by simulator + decision_var dispatcher)
global use_improved
use_improved = logical(cfg.use_improved);

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
    results = run_benchmark_compare(cfg); %#ok<NASGU>
    
    % Generate smart filename: benchmark_results_<solver>_<variant>_<strategies>.mat
    solvers = string(cfg.bench_solvers);
    strategies = string(cfg.bench_strategies);
    solver_str = char(solvers(1));  % typically just one solver in benchmark
    strategy_str = strjoin(strategies, '-');  % 'base-improved'
    
    % Add variant info if nsga2
    if strcmp(solver_str, 'nsga2')
        variant_str = char(cfg.nsga2_variant);
        benchmark_fname = sprintf('benchmark_results_%s_%s_%s.mat', ...
            solver_str, variant_str, strategy_str);
    else
        benchmark_fname = sprintf('benchmark_results_%s_%s.mat', ...
            solver_str, strategy_str);
    end
    
    save(benchmark_fname, 'results');
    fprintf('Saved: %s\n', benchmark_fname);
else
    OUT = run_single(cfg); %#ok<NASGU>
    
    % Generate smart filename based on optimizer + variant + strategy
    segment = extract_segment_name(cfg.route_file);
    result_fname = sprintf('result_%s_%s_%s%s.mat', ...
        OUT.optimizer, OUT.nsga2_variant, OUT.strategy, segment);
    
    save(result_fname, 'OUT');
    fprintf('Saved: %s\n', result_fname);
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

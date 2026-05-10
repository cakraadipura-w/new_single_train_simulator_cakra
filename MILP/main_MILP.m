%% main_MILP.m - MILP train control prototype
% Run this file directly from the MATLAB editor.
%
% Notes:
% - Requires YALMIP and a MILP-capable backend solver.
% - This script auto-detects a local YALMIP checkout and can reuse the
%   MOSEK MATLAB binary bundled under a CVX installation.

clearvars -except route_direction selected_segments T_target_override dx_nominal pwl_segment_count milp_solver_name show_figures verbose_level
clc; close all;

%% ---------------------- User settings ----------------------
if ~exist('route_direction', 'var') || isempty(route_direction)
    route_direction = 'up';       % 'up' | 'down'
end
if ~exist('selected_segments', 'var') || isempty(selected_segments)
    selected_segments = {'IS01'}; % {'IS01'}, {'IS01','IS02'}, or 'all'
end
if ~exist('T_target_override', 'var') || isempty(T_target_override)
    T_target_override = [];       % [] -> use scheduled time from catalog
end
if ~exist('dx_nominal', 'var') || isempty(dx_nominal)
    dx_nominal = 50;              % nominal space step (m), start coarse for faster MILP runs
end
if ~exist('pwl_segment_count', 'var') || isempty(pwl_segment_count)
    pwl_segment_count = 8;        % number of PWL intervals, increase later for finer fidelity
end
if ~exist('milp_solver_name', 'var') || isempty(milp_solver_name)
    milp_solver_name = 'sdpt3';   % 'mosek', 'gurobi', 'sdpt3' or another YALMIP MILP solver
end
if ~exist('show_figures', 'var') || isempty(show_figures)
    show_figures = true;
end
if ~exist('verbose_level', 'var') || isempty(verbose_level)
    verbose_level = 1;
end

%% ---------------------- Setup ----------------------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(genpath(project_root));

ensure_yalmip_ready(project_root, script_dir);
milp_solver_name = prepare_milp_solver(milp_solver_name, project_root, script_dir);
solver_tag = normalize_solver_token(milp_solver_name);

segment_catalog = build_segment_catalog(route_direction);
selected_catalog = resolve_segment_selection(selected_segments, segment_catalog);

results_root = fullfile(script_dir, 'results');
if ~exist(results_root, 'dir')
    mkdir(results_root);
end

fprintf('Loading rolling stock...\n');
rollingstock_file = fullfile(project_root, 'rollingstocks', 'rollingstock_Guangzhou_L7.m');
run(rollingstock_file);

params = struct();
params.gravity = 9.81;
params.M_eff = Mass * (1 + lambda) * 1000;      % kg
params.Davis = Davis;
params.max_traction_N = max_traction * 1000;    % N
params.max_accel_trac = max_accel_trac;
params.max_accel_brake = max_accel_brake;
params.Max_tractive_power = Max_tractive_power; % W
params.min_avg_speed = 0.10;                    % m/s, reciprocal guard

results = repmat(init_result_struct(), 0, 1);

%% ---------------------- Solve per segment ----------------------
for seg_idx = 1:numel(selected_catalog)
    seg = selected_catalog(seg_idx);

    if isempty(T_target_override)
        T_target = seg.T_sched;
    else
        T_target = T_target_override;
    end

    route_file = resolve_project_route_path(project_root, seg.file, route_direction);
    out_dir = fullfile(results_root, seg.name);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    fprintf('\n============================================================\n');
    fprintf('MILP optimization   | %s\n', seg.name);
    fprintf('Route file          | %s\n', route_file);
    fprintf('Target time         | %.2f s\n', T_target);
    fprintf('Results folder      | %s\n', out_dir);
    fprintf('============================================================\n');

    segment_timer = tic;
    route_data = load(route_file);
    problem = build_segment_problem(route_data, params, T_target, dx_nominal, pwl_segment_count);
    solution = solve_milp_segment(problem, params, milp_solver_name, verbose_level);
    processing_time_s = toc(segment_timer);

    result = build_segment_result(seg, route_file, out_dir, problem, solution, solver_tag, processing_time_s);
    save_segment_outputs(result, show_figures);
    results(end+1, 1) = result; %#ok<AGROW>
end

batch_file = fullfile(results_root, sprintf('MILP_batch_results_%s.mat', solver_tag));
save(batch_file, 'results');
fprintf('\nSaved batch results: %s\n', batch_file);

%% ---------------------- Local helpers ----------------------

function catalog = build_segment_catalog(route_direction)
    catalog = get_guangzhou_line7_catalog(route_direction);
end

function selected_catalog = resolve_segment_selection(selected_segments, catalog)
    if ischar(selected_segments) || (isstring(selected_segments) && isscalar(selected_segments))
        selected_segments = cellstr(string(selected_segments));
    elseif isstring(selected_segments)
        selected_segments = cellstr(selected_segments(:));
    elseif ~iscell(selected_segments)
        error('selected_segments must be ''all'', a string, or a cell array of strings.');
    end

    if isscalar(selected_segments) && strcmpi(string(selected_segments{1}), 'all')
        selected_catalog = catalog;
        return;
    end

    selected_catalog = repmat(catalog(1), 0, 1);
    keep_mask = false(1, numel(catalog));
    names = {catalog.name};
    files = {catalog.file};

    for idx = 1:numel(selected_segments)
        token = char(string(selected_segments{idx}));
        match_idx = find(strcmpi(names, token) | strcmpi(files, token), 1);
        if isempty(match_idx)
            error('Unknown segment selection: %s', token);
        end
        if ~keep_mask(match_idx)
            selected_catalog(end+1, 1) = catalog(match_idx); %#ok<AGROW>
            keep_mask(match_idx) = true;
        end
    end
end

function ensure_yalmip_ready(project_root, script_dir)
    if exist('sdpvar', 'file') == 2 && exist('optimize', 'file') == 2
        fprintf('YALMIP detected on MATLAB path.\n');
        return;
    end

    candidate_roots = build_yalmip_candidate_roots(project_root, script_dir);
    yalmip_root = '';

    for idx = 1:numel(candidate_roots)
        yalmip_root = find_yalmip_root(candidate_roots{idx});
        if ~isempty(yalmip_root)
            break;
        end
    end

    if isempty(yalmip_root)
        error(['YALMIP was not found on the MATLAB path or in common install folders.' newline ...
            'Install steps:' newline ...
            '1. Download or clone YALMIP (https://github.com/yalmip/YALMIP).' newline ...
            '2. Put it in a folder such as <project>/third_party/YALMIP.' newline ...
            '3. Then rerun main_MILP.m']);
    end

    addpath(genpath(yalmip_root));

    if exist('sdpvar', 'file') ~= 2 || exist('optimize', 'file') ~= 2
        error('YALMIP path setup ran but sdpvar/optimize are still unavailable.');
    end

    fprintf('YALMIP ready at: %s\n', yalmip_root);
end

function solver_name = prepare_milp_solver(solver_name, project_root, script_dir)
    solver_name = strtrim(char(string(solver_name)));
    if isempty(solver_name)
        error('milp_solver_name must not be empty.');
    end

    requested_solver = lower(solver_name);
    unsupported_milp_solvers = {'sdpt3', 'sedumi', 'ecos', 'ecos_bb', 'scs', 'osqp', 'quadprog', 'linprog'};
    if any(strcmp(requested_solver, unsupported_milp_solvers))
        fallback_solver = pick_available_milp_solver(project_root, script_dir);
        fprintf(['MILP solver selection | %s is not a MILP solver for this model.' newline ...
            'Falling back to %s.\n'], solver_name, fallback_solver);
        solver_name = fallback_solver;
        requested_solver = lower(solver_name);
    end

    switch requested_solver
        case 'mosek'
            ensure_solver_binary('mosekopt', build_mosek_candidate_dirs(project_root, script_dir), ...
                ['MOSEK MATLAB binary was not found.' newline ...
                'Expected a mosekopt.mexw64 file, for example under:' newline ...
                'C:\Program Files\MATLAB\cvx\mosek\w64']);
        case 'gurobi'
            ensure_solver_binary('gurobi', {}, 'Gurobi MATLAB binary was not found on the MATLAB path.');
        case {'cplex', 'cplexmilp'}
            if exist('cplexmilp', 'file') ~= 2 && exist('cplexlp', 'file') ~= 2
                error('CPLEX MATLAB binaries were not found on the MATLAB path.');
            end
        otherwise
            fallback_solver = pick_available_milp_solver(project_root, script_dir, requested_solver);
            if ~strcmpi(fallback_solver, solver_name)
                fprintf(['MILP solver selection | requested solver %s is not ready for MILP in this MATLAB session.' newline ...
                    'Falling back to %s.\n'], solver_name, fallback_solver);
                solver_name = fallback_solver;
            else
                fprintf('MILP solver selection | relying on YALMIP solver name: %s\n', solver_name);
                return;
            end
    end

    fprintf('MILP solver selection | %s\n', solver_name);
end

function solver_name = pick_available_milp_solver(project_root, script_dir, preferred_solver)
    if nargin < 3
        preferred_solver = '';
    end

    candidates = {};
    if ~isempty(preferred_solver)
        candidates{end+1, 1} = lower(strtrim(char(string(preferred_solver)))); %#ok<AGROW>
    end
    candidates = [candidates; {'mosek'; 'gurobi'; 'cplex'}];

    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    ordered = cell(0, 1);
    for idx = 1:numel(candidates)
        name = candidates{idx};
        if isempty(name) || isKey(seen, name)
            continue;
        end
        seen(name) = true;
        ordered{end+1, 1} = name; %#ok<AGROW>
    end

    for idx = 1:numel(ordered)
        switch ordered{idx}
            case 'mosek'
                if solver_binary_exists('mosekopt', build_mosek_candidate_dirs(project_root, script_dir), true)
                    solver_name = 'mosek';
                    return;
                end
            case 'gurobi'
                if solver_binary_exists('gurobi', {}, false)
                    solver_name = 'gurobi';
                    return;
                end
            case 'cplex'
                if exist('cplexmilp', 'file') == 2 || exist('cplexlp', 'file') == 2
                    solver_name = 'cplex';
                    return;
                end
        end
    end

    error(['No MILP-capable YALMIP backend is available.' newline ...
        'Install or enable one of: MOSEK, Gurobi, or CPLEX.']);
end

function ensure_solver_binary(binary_name, candidate_dirs, error_message)
    if solver_binary_exists(binary_name, candidate_dirs, true)
        return;
    end

    error(error_message);
end

function tf = solver_binary_exists(binary_name, candidate_dirs, addCandidatesToPath)
    if nargin < 3
        addCandidatesToPath = true;
    end

    if exist(binary_name, 'file') == 2 || exist(binary_name, 'file') == 3
        tf = true;
        return;
    end

    for idx = 1:numel(candidate_dirs)
        candidate_dir = candidate_dirs{idx};
        if isempty(candidate_dir) || ~exist(candidate_dir, 'dir')
            continue;
        end
        if addCandidatesToPath
            addpath(candidate_dir);
        end
        if exist(binary_name, 'file') == 2 || exist(binary_name, 'file') == 3
            tf = true;
            return;
        end
    end

    tf = false;
end

function candidate_roots = build_yalmip_candidate_roots(project_root, script_dir)
    user_home = getenv('USERPROFILE');
    if isempty(user_home)
        user_home = char(java.lang.System.getProperty('user.home'));
    end

    candidate_roots = { ...
        fullfile(project_root, 'third_party', 'YALMIP'), ...
        fullfile(project_root, 'YALMIP'), ...
        fullfile(script_dir, 'YALMIP'), ...
        fullfile(user_home, 'Documents', 'MATLAB', 'YALMIP'), ...
        fullfile(user_home, 'MATLAB', 'YALMIP'), ...
        fullfile(user_home, 'Downloads', 'YALMIP'), ...
        fullfile(user_home, 'Downloads', 'yalmip-master'), ...
        fullfile(user_home, 'Desktop', 'YALMIP')};

    candidate_roots = unique(candidate_roots(:), 'stable');
end

function yalmip_root = find_yalmip_root(candidate)
    yalmip_root = '';

    if isempty(candidate) || ~exist(candidate, 'dir')
        return;
    end

    if exist(fullfile(candidate, '@sdpvar', 'sdpvar.m'), 'file') == 2
        yalmip_root = candidate;
        return;
    end

    entries = dir(candidate);
    for idx = 1:numel(entries)
        if ~entries(idx).isdir
            continue;
        end
        name = entries(idx).name;
        if strcmp(name, '.') || strcmp(name, '..')
            continue;
        end
        child = fullfile(candidate, name);
        if exist(fullfile(child, '@sdpvar', 'sdpvar.m'), 'file') == 2
            yalmip_root = child;
            return;
        end
    end
end

function candidate_dirs = build_mosek_candidate_dirs(project_root, script_dir)
    candidate_dirs = {};

    cvx_setup_path = which('cvx_setup');
    if ~isempty(cvx_setup_path)
        candidate_dirs{end+1, 1} = fullfile(fileparts(cvx_setup_path), 'mosek', 'w64'); %#ok<AGROW>
    end

    candidate_dirs{end+1, 1} = fullfile(project_root, 'cvx', 'mosek', 'w64'); %#ok<AGROW>
    candidate_dirs{end+1, 1} = fullfile(script_dir, 'cvx', 'mosek', 'w64'); %#ok<AGROW>

    env_names = {'ProgramFiles', 'ProgramW6432', 'ProgramFiles(x86)'};
    for idx = 1:numel(env_names)
        base_dir = getenv(env_names{idx});
        if isempty(base_dir)
            continue;
        end
        candidate_dirs{end+1, 1} = fullfile(base_dir, 'MATLAB', 'cvx', 'mosek', 'w64'); %#ok<AGROW>

        matlab_root = fullfile(base_dir, 'MATLAB');
        if ~exist(matlab_root, 'dir')
            continue;
        end
        entries = dir(matlab_root);
        for entry_idx = 1:numel(entries)
            if ~entries(entry_idx).isdir
                continue;
            end
            name = entries(entry_idx).name;
            if strcmp(name, '.') || strcmp(name, '..')
                continue;
            end
            candidate_dirs{end+1, 1} = fullfile(matlab_root, name, 'cvx', 'mosek', 'w64'); %#ok<AGROW>
        end
    end

    candidate_dirs = unique(candidate_dirs, 'stable');
end

function solver_tag = normalize_solver_token(solver_name)
    solver_tag = strtrim(char(string(solver_name)));
    if isempty(solver_tag)
        solver_tag = 'default';
        return;
    end

    solver_tag = lower(regexprep(solver_tag, '[^A-Za-z0-9]+', '_'));
    solver_tag = regexprep(solver_tag, '^_+|_+$', '');
    if isempty(solver_tag)
        solver_tag = 'default';
    end
end

function problem = build_segment_problem(route_data, params, T_target, dx_nominal, pwl_segment_count)
    if ~isfield(route_data, 'vel_profile')
        error('Route file is missing vel_profile.');
    end

    vel_lim_raw = route_data.vel_profile;
    if isfield(route_data, 'gradient')
        grad_raw = route_data.gradient;
    else
        grad_raw = [vel_lim_raw(:,1), zeros(size(vel_lim_raw,1), 1)];
    end

    s_total = vel_lim_raw(end, 1) * 1000;
    s_nodes = (0:dx_nominal:s_total)';
    if abs(s_nodes(end) - s_total) > 1e-9
        s_nodes = [s_nodes; s_total]; %#ok<AGROW>
    end

    ds = diff(s_nodes);
    s_mid = s_nodes(1:end-1) + 0.5 * ds;

    v_lim_nodes = interp1(vel_lim_raw(:,1) * 1000, vel_lim_raw(:,2) / 3.6, ...
        s_nodes, 'previous', 'extrap');
    g_permil_mid = interp1(grad_raw(:,1) * 1000, grad_raw(:,2), s_mid, 'linear', 0);

    v_max = max(v_lim_nodes) + 2;
    v_break = linspace(0, v_max, pwl_segment_count + 1)';
    v_break_sq = v_break .^ 2;
    v_break_inv = 1 ./ max(v_break, params.min_avg_speed);

    problem = struct();
    problem.T_target = T_target;
    problem.s_nodes = s_nodes;
    problem.ds = ds;
    problem.s_mid = s_mid;
    problem.N = numel(s_nodes);
    problem.M = numel(ds);
    problem.J = pwl_segment_count;
    problem.v_lim_nodes = v_lim_nodes(:);
    problem.G_seg = params.M_eff * params.gravity * (g_permil_mid(:) / 1000);
    problem.v_break = v_break;
    problem.v_break_sq = v_break_sq;
    problem.v_break_inv = v_break_inv;
end

function solution = solve_milp_segment(problem, params, solver_name, verbose_level)
    N = problem.N;
    M = problem.M;
    ds = problem.ds;
    v_lim_nodes = problem.v_lim_nodes;
    G_seg = problem.G_seg;
    v_break = problem.v_break;
    v_break_sq = problem.v_break_sq;
    v_break_inv = problem.v_break_inv;

    davis_A = params.Davis(1) * 1000;
    davis_B = params.Davis(2) * 1000;
    davis_C = params.Davis(3) * 1000;

    v = sdpvar(N, 1);
    speed_sq = sdpvar(N, 1);
    lambda = sdpvar(N, numel(v_break), 'full');

    v_avg = sdpvar(M, 1);
    alpha = sdpvar(M, 1);
    mu = sdpvar(M, numel(v_break), 'full');

    f_trac = sdpvar(M, 1);
    t_seg = sdpvar(M, 1);

    Constraints = [];

    for i = 1:N
        Constraints = [Constraints, ...
            sum(lambda(i,:)) == 1, ...
            lambda(i,:) >= 0, ...
            v(i) == lambda(i,:) * v_break, ...
            speed_sq(i) == lambda(i,:) * v_break_sq, ...
            sos2(lambda(i,:))];
    end

    for i = 1:M
        Constraints = [Constraints, ...
            v_avg(i) == 0.5 * (v(i) + v(i+1)), ...
            sum(mu(i,:)) == 1, ...
            mu(i,:) >= 0, ...
            v_avg(i) == mu(i,:) * v_break, ...
            alpha(i) == mu(i,:) * v_break_inv, ...
            sos2(mu(i,:))];
    end

    Constraints = [Constraints, ...
        v(1) == 0, ...
        v(N) == 0, ...
        v >= 0, ...
        v <= v_lim_nodes, ...
        f_trac >= 0, ...
        v_avg >= params.min_avg_speed, ...
        t_seg == ds .* alpha, ...
        sum(t_seg) <= problem.T_target];

    for i = 1:M
        beta_avg = 0.5 * (speed_sq(i) + speed_sq(i+1));
        acc_i = (speed_sq(i+1) - speed_sq(i)) / (2 * ds(i));
        resist_i = davis_A + davis_B * v_avg(i) + davis_C * beta_avg;

        Constraints = [Constraints, ...
            f_trac(i) >= params.M_eff * acc_i + G_seg(i) + resist_i, ...
            f_trac(i) <= params.max_traction_N, ...
            f_trac(i) <= params.Max_tractive_power * alpha(i), ...
            acc_i <= params.max_accel_trac, ...
            acc_i >= -params.max_accel_brake];
    end

    Objective = sum(f_trac .* ds);

    ops = sdpsettings('solver', solver_name, 'verbose', verbose_level);
    solve_timer = tic;
    diagnosis = optimize(Constraints, Objective, ops);
    solve_wall_time_s = toc(solve_timer);
    if diagnosis.problem ~= 0
        error('MILP solver failed: %s', diagnosis.info);
    end

    solution = struct();
    solution.status = 'Solved';
    solution.diagnosis_info = diagnosis.info;
    solution.optval = value(Objective);
    solution.v = full(value(v));
    solution.speed_sq = full(value(speed_sq));
    solution.v_avg = full(value(v_avg));
    solution.alpha = full(value(alpha));
    solution.f_trac = full(value(f_trac));
    solution.dt = full(value(t_seg));
    solution.T_actual = sum(solution.dt);
    solution.E_kWh = solution.optval / 3.6e6;
    solution.solve_time_s = extract_diagnostic_seconds(diagnosis, 'solvertime', solve_wall_time_s);
end

function result = build_segment_result(seg, route_file, out_dir, problem, solution, solver_tag, processing_time_s)
    result = init_result_struct();
    result.segment = seg.name;
    result.solver = solver_tag;
    result.route_file = route_file;
    result.T_target = problem.T_target;
    result.status = solution.status;
    result.diagnosis_info = solution.diagnosis_info;
    result.optval = solution.optval;
    result.s_nodes = problem.s_nodes;
    result.s_mid = problem.s_mid;
    result.v_mps = solution.v;
    result.v_kmph = solution.v * 3.6;
    result.v_lim_kmph = problem.v_lim_nodes * 3.6;
    result.f_trac_N = solution.f_trac;
    result.f_trac_kN = solution.f_trac / 1000;
    result.dt = solution.dt;
    result.T_actual = solution.T_actual;
    result.E_kWh = solution.E_kWh;
    result.solve_time_s = solution.solve_time_s;
    result.processing_time_s = processing_time_s;
    result.output_dir = out_dir;
    result.result_mat = fullfile(out_dir, sprintf('MILP_%s_%s_result.mat', seg.name, solver_tag));
    result.summary_txt = fullfile(out_dir, sprintf('MILP_%s_%s_summary.txt', seg.name, solver_tag));
    result.plot_png = fullfile(out_dir, sprintf('MILP_%s_%s_profiles.png', seg.name, solver_tag));
end

function save_segment_outputs(result, show_figures)
    save(result.result_mat, 'result');

    fid = fopen(result.summary_txt, 'w');
    if fid ~= -1
        fprintf(fid, 'Segment              : %s\n', result.segment);
        fprintf(fid, 'Solver               : %s\n', result.solver);
        fprintf(fid, 'Route file           : %s\n', result.route_file);
        fprintf(fid, 'Status               : %s\n', result.status);
        fprintf(fid, 'Solver info          : %s\n', result.diagnosis_info);
        fprintf(fid, 'Target time          : %s\n', format_time_seconds_minutes(result.T_target));
        fprintf(fid, 'Actual time          : %s\n', format_time_seconds_minutes(result.T_actual));
        fprintf(fid, 'Energy (kWh)         : %.6f\n', result.E_kWh);
        fprintf(fid, 'Solve time           : %s\n', format_time_seconds_minutes(result.solve_time_s));
        fprintf(fid, 'Processing time      : %s\n', format_time_seconds_minutes(result.processing_time_s));
        fprintf(fid, 'Result MAT           : %s\n', result.result_mat);
        fprintf(fid, 'Profile plot         : %s\n', result.plot_png);
        fclose(fid);
    end

    fig = figure('Visible', ternary(show_figures, 'on', 'off'), 'Color', 'w', ...
        'Name', sprintf('MILP | %s', result.segment), 'NumberTitle', 'off');

    subplot(2,1,1);
    plot(result.s_nodes / 1000, result.v_kmph, 'b-', 'LineWidth', 2); hold on;
    stairs(result.s_nodes / 1000, result.v_lim_kmph, 'k--', 'LineWidth', 1.3);
    xlabel('Distance (km)');
    ylabel('Speed (km/h)');
    title(sprintf('%s | %s | status=%s | T=%.2f s | E=%.3f kWh', ...
        result.segment, upper(result.solver), result.status, result.T_actual, result.E_kWh));
    legend('MILP optimal speed', 'Speed limit', 'Location', 'best');
    grid on;

    subplot(2,1,2);
    plot(result.s_mid / 1000, result.f_trac_kN, 'r-', 'LineWidth', 1.5);
    xlabel('Distance (km)');
    ylabel('Tractive force (kN)');
    title('Traction force profile');
    grid on;

    saveas(fig, result.plot_png);
    if ~show_figures
        close(fig);
    end

    fprintf('Status: %s | T_actual=%s | E=%.3f kWh | solve=%s | process=%s\n', ...
        result.status, format_time_seconds_minutes(result.T_actual), result.E_kWh, ...
        format_time_seconds_minutes(result.solve_time_s), format_time_seconds_minutes(result.processing_time_s));
    fprintf('Saved result  : %s\n', result.result_mat);
    fprintf('Saved summary : %s\n', result.summary_txt);
    fprintf('Saved plot    : %s\n', result.plot_png);
end

function result = init_result_struct()
    result = struct( ...
        'segment', '', ...
        'solver', '', ...
        'route_file', '', ...
        'T_target', NaN, ...
        'status', '', ...
        'diagnosis_info', '', ...
        'optval', NaN, ...
        's_nodes', zeros(0,1), ...
        's_mid', zeros(0,1), ...
        'v_mps', zeros(0,1), ...
        'v_kmph', zeros(0,1), ...
        'v_lim_kmph', zeros(0,1), ...
        'f_trac_N', zeros(0,1), ...
        'f_trac_kN', zeros(0,1), ...
        'dt', zeros(0,1), ...
        'T_actual', NaN, ...
        'E_kWh', NaN, ...
        'solve_time_s', NaN, ...
        'processing_time_s', NaN, ...
        'output_dir', '', ...
        'result_mat', '', ...
        'summary_txt', '', ...
        'plot_png', '');
end

function seconds = extract_diagnostic_seconds(diagnosis, field_name, fallback_value)
    seconds = fallback_value;
    if isstruct(diagnosis) && isfield(diagnosis, field_name)
        candidate = diagnosis.(field_name);
        if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
            seconds = double(candidate);
        end
    end
end

function text = format_time_seconds_minutes(seconds)
    text = sprintf('%.6f s (%.6f min)', seconds, seconds / 60);
end

function out = ternary(cond, true_value, false_value)
    if cond
        out = true_value;
    else
        out = false_value;
    end
end
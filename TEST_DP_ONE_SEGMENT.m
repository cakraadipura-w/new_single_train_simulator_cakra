%% =========================================================================
% TEST_DP_ONE_SEGMENT_V3.m
% Perbaikan unit: Davis dalam kN (v dalam m/s), biaya DP dalam kWh
% =========================================================================
clear; clc; close all;

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(project_root));

ALL_SEGMENTS = struct( ...
    'name',    {'IS01','IS02','IS03','IS04','IS05','IS06','IS07','IS08'}, ...
    'file',    { ...
        'Guangzhou_Line7_IS01_0.000-1.120km.mat', ...
        'Guangzhou_Line7_IS02_1.120-3.028km.mat', ...
        'Guangzhou_Line7_IS03_3.028-5.200km.mat', ...
        'Guangzhou_Line7_IS04_5.200-6.842km.mat', ...
        'Guangzhou_Line7_IS05_6.842-8.958km.mat', ...
        'Guangzhou_Line7_IS06_8.958-11.323km.mat', ...
        'Guangzhou_Line7_IS07_11.323-13.729km.mat', ...
        'Guangzhou_Line7_IS08_13.729-17.507km.mat'}, ...
    'T_sched', {129, 169, 184, 177, 185, 219, 214, 329});

selected_segments = 'all';  % bisa juga 'all', 1..8, [1 3 5], {'IS01','IS03'}
segments_to_run = resolve_selected_segments(selected_segments, ALL_SEGMENTS);

target_T_offsets = 0:4;  % otomatis jadi [T_sched, T_sched+1, ...]
target_T_values = [];    % isi manual kalau mau override, contoh: [129, 130, 131]
time_tolerance_s = 0.5;        % toleransi kecocokan akhir terhadap target
control_refinement_steps = [0.15, 0.05, 0.025];
lambda_bounds = [0.005, 0.2];
lambda_search_iters = 8;
lambda_seed_count = 2;
auto_adjust_dv = true;
min_control_step_small_dx = 0.05;
max_lambda_expand_steps = 10;
min_lambda_value = 1e-4;
max_lambda_value = 10;
summary_output_dir = fullfile(project_root, 'dp_target_time_results');
if exist(summary_output_dir, 'dir') ~= 7
    mkdir(summary_output_dir);
end

%% 1. Load Rolling Stock
fprintf('Loading rolling stock...\n');
run(fullfile(project_root, 'rollingstocks', 'rollingstock_Guangzhou_L7.m'));
inertial_mass_kg = Mass * (1 + lambda) * 1000;
overall_batch_timer = tic;
segment_runs = repmat(struct( ...
    'segment', struct('name', '', 'file', '', 'T_sched', NaN), ...
    'target_T_values', [], ...
    'target_results', [], ...
    'total_search_time_s', NaN, ...
    'summary_file', '', ...
    'results_file', '', ...
    'run_meta', struct()), 1, numel(segments_to_run));

for seg_run_idx = 1:numel(segments_to_run)
    active_segment = segments_to_run(seg_run_idx);
    segment_target_T_values = resolve_target_values(target_T_values, active_segment.T_sched, target_T_offsets);

    [target_results, total_search_time_s, run_meta] = run_segment_target_search( ...
        active_segment, segment_target_T_values, ...
        time_tolerance_s, control_refinement_steps, lambda_bounds, ...
        lambda_search_iters, lambda_seed_count, auto_adjust_dv, ...
        min_control_step_small_dx, max_lambda_expand_steps, ...
        min_lambda_value, max_lambda_value, project_root, ...
        Mass, gravity, Davis, inertial_mass_kg, ...
        Max_tractive_power, V1_traction, V2_traction, ...
        Max_brake_power, V1_brake, V2_brake);

    segment_summary_text = build_segment_summary_text( ...
        active_segment, segment_target_T_values, target_results, total_search_time_s, run_meta);
    fprintf('\n%s', segment_summary_text);

    route_base_name = erase(active_segment.file, '.mat');
    segment_summary_file = fullfile(summary_output_dir, sprintf('%s_summary.txt', route_base_name));
    segment_results_file = fullfile(summary_output_dir, sprintf('%s_target_time_results.mat', route_base_name));
    write_text_file(segment_summary_file, segment_summary_text);

    active_segment_result = active_segment;
    target_values_saved = segment_target_T_values;
    total_search_time_saved = total_search_time_s;
    run_meta_saved = run_meta;
    save(segment_results_file, 'active_segment_result', 'target_values_saved', ...
        'target_results', 'total_search_time_saved', 'run_meta_saved');

    segment_runs(seg_run_idx) = struct( ...
        'segment', active_segment, ...
        'target_T_values', segment_target_T_values, ...
        'target_results', target_results, ...
        'total_search_time_s', total_search_time_s, ...
        'summary_file', segment_summary_file, ...
        'results_file', segment_results_file, ...
        'run_meta', run_meta);
end

overall_batch_time_s = toc(overall_batch_timer);
aggregate_summary_text = build_overall_summary_text(segment_runs, overall_batch_time_s, summary_output_dir);
fprintf('\n%s', aggregate_summary_text);

aggregate_summary_file = fullfile(summary_output_dir, 'all_segments_target_time_summary.txt');
aggregate_results_file = fullfile(summary_output_dir, 'all_segments_target_time_results.mat');
write_text_file(aggregate_summary_file, aggregate_summary_text);
save(aggregate_results_file, 'segment_runs', 'overall_batch_time_s', 'segments_to_run', ...
    'target_T_offsets', 'target_T_values');

fprintf('Summary files saved in: %s\n', summary_output_dir);

%% 3. Fungsi Fisika
function F = traction_force(v, Pmax, V1, V2)
    % Returns traction force in N
    if v <= V1
        F = Pmax / V1;
    elseif v <= V2
        F = Pmax / v;
    else
        F = Pmax * V2 / v^2;
    end
end

function F = brake_force(v, Pmax, V1, V2)
    % Returns braking force magnitude in N
    if v <= V1
        F = Pmax / V1;
    elseif v <= V2
        F = Pmax / v;
    else
        F = Pmax * V2 / v^2;
    end
end

%% 4. DP Solver
% Cost function: J = E_kWh + lambda * T_s
% lambda in kWh/s — typical range 1e-4 to 1e-2 for this segment
function [T_total, E_total_kWh, v_opt, t_opt, x_opt] = solve_dp(lambda, ...
    vel_lim, grad, x_grid, v_grid, control_levels, dx, a_max, a_min, epsilon_v, ...
    Mass, gravity, Davis, inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

    N_stages = length(x_grid);
    Nv       = length(v_grid);
    Nc       = length(control_levels);

    J      = inf(Nv, N_stages);
    policy = nan(Nv, N_stages-1);

    % Terminal: all states with v <= epsilon_v cost 0
    for iv = 1:Nv
        if v_grid(iv) <= epsilon_v
            J(iv, end) = 0;
        end
    end

    % Backward DP
    for k = N_stages-1 : -1 : 1
        x_cur       = x_grid(k);
        v_lim       = interp1(vel_lim(:,1), vel_lim(:,2), x_cur, 'previous', vel_lim(end,2));
        grad_permil = interp1(grad(:,1), grad(:,2), x_cur, 'previous', 0);
        theta       = atan(grad_permil / 1000);
        F_grade     = Mass * 1000 * gravity * sin(theta);   % N

        for iv = 1:Nv
            v = v_grid(iv);
            if v > v_lim + 1e-6
                continue;
            end

            % Davis resistance in N: Davis coefficients in kN, v in m/s
            R_davis = (Davis(1) + Davis(2)*v + Davis(3)*v^2) * 1000;  % N

            for ic = 1:Nc
                c = control_levels(ic);
                if c >= 0
                    F_max     = traction_force(v, Max_tractive_power, V1_traction, V2_traction);
                    F_applied = c * F_max;
                else
                    F_max     = brake_force(v, Max_brake_power, V1_brake, V2_brake);
                    F_applied = c * F_max;   % negative
                end

                F_net = F_applied - R_davis - F_grade;
                a     = F_net / inertial_mass_kg;
                if a > a_max || a < a_min
                    continue;
                end

                v_sq = v^2 + 2*a*dx;
                if v_sq < 0
                    continue;
                end
                v_next = sqrt(v_sq);

                dt    = 2*dx / (v + max(v_next, 0.01));
                % Energy in kWh (only traction counts)
                dE_kWh = max(0, F_applied) * dx / 3.6e6;

                if k == N_stages-1
                    if v_next <= epsilon_v
                        J_cand = dE_kWh + lambda*dt;
                        if J_cand < J(iv, k)
                            J(iv, k)      = J_cand;
                            policy(iv, k) = c;
                        end
                    end
                else
                    J_next = interpolate_cost_to_go(v_next, v_grid, J(:, k+1));
                    if isinf(J_next)
                        continue;
                    end
                    J_cand = dE_kWh + lambda*dt + J_next;
                    if J_cand < J(iv, k)
                        J(iv, k)      = J_cand;
                        policy(iv, k) = c;
                    end
                end
            end

            if k == N_stages-1 && v > epsilon_v
                a_stop = -v^2 / (2*dx);
                if a_stop >= a_min - 1e-9
                    F_brake_cap = brake_force(v, Max_brake_power, V1_brake, V2_brake);
                    F_applied_stop = inertial_mass_kg * a_stop + R_davis + F_grade;
                    if F_applied_stop <= 0 && F_applied_stop >= -F_brake_cap - 1e-9
                        c_stop = max(-1, F_applied_stop / max(F_brake_cap, 1));
                        dt_stop = 2*dx / max(v, 0.01);
                        J_cand = lambda * dt_stop;
                        if J_cand < J(iv, k)
                            J(iv, k)      = J_cand;
                            policy(iv, k) = c_stop;
                        end
                    end
                end
            end
        end
    end

    % Diagnosis: check if start state is reachable
    [~, iv0] = min(abs(v_grid));
    if isinf(J(iv0, 1))
        warning('solve_dp: J(v=0, stage=1) = inf. No feasible path found (lambda=%.4g).', lambda);
    end

    % Forward propagation
    iv    = iv0;
    v_opt = zeros(1, N_stages);
    t_opt = zeros(1, N_stages);
    E_opt = zeros(1, N_stages);
    x_opt = x_grid;

    for k = 1:N_stages-1
        c = policy(iv, k);
        if isnan(c)
            iv_fallback = nearest_feasible_policy_bin(policy(:, k), v_opt(k), v_grid);
            if iv_fallback == 0
                warning('solve_dp: no feasible control at stage %d (iv=%d). Stopping.', k, iv);
                break;
            end
            c = policy(iv_fallback, k);
        end
        v     = v_opt(k);
        x_cur = x_grid(k);

        grad_permil = interp1(grad(:,1), grad(:,2), x_cur, 'previous', 0);
        theta   = atan(grad_permil/1000);
        F_grade = Mass*1000*gravity*sin(theta);
        R_davis = (Davis(1) + Davis(2)*v + Davis(3)*v^2) * 1000;  % N

        if c >= 0
            F_max     = traction_force(v, Max_tractive_power, V1_traction, V2_traction);
            F_applied = c * F_max;
        else
            F_max     = brake_force(v, Max_brake_power, V1_brake, V2_brake);
            F_applied = c * F_max;
        end

        a      = (F_applied - R_davis - F_grade) / inertial_mass_kg;
        v_sq   = v^2 + 2*a*dx;
        v_next = sqrt(max(v_sq, 0));

        if k == N_stages-1
            v_next = min(v_next, epsilon_v);
        end
        [~, iv] = min(abs(v_grid - v_next));
        v_opt(k+1) = v_next;
        dt         = 2*dx / (v + max(v_next, 0.01));
        t_opt(k+1) = t_opt(k) + dt;
        E_opt(k+1) = E_opt(k) + max(0, F_applied)*dx / 3.6e6;  % kWh
    end

    T_total      = t_opt(end);
    E_total_kWh  = E_opt(end);
end

function active_segments = resolve_selected_segments(selected_segments, all_segments)
    if isnumeric(selected_segments)
        run_idx = selected_segments(:)';
    else
        selected_tokens = cellstr(string(selected_segments));
        if numel(selected_tokens) == 1 && strcmpi(selected_tokens{1}, 'all')
            run_idx = 1:numel(all_segments);
        else
            run_idx = zeros(1, numel(selected_tokens));
            for idx = 1:numel(selected_tokens)
                seg_idx = find(strcmpi({all_segments.name}, selected_tokens{idx}), 1, 'first');
                assert(~isempty(seg_idx), ...
                    'Segment tidak valid: %s. Gunakan IS01..IS08 atau ''all''.', selected_tokens{idx});
                run_idx(idx) = seg_idx;
            end
        end
    end

    assert(~isempty(run_idx), 'Pilih minimal satu segment untuk dijalankan.');
    run_idx = unique(run_idx, 'stable');
    assert(all(run_idx >= 1 & run_idx <= numel(all_segments)), ...
        'Index segment harus berada pada rentang 1..%d.', numel(all_segments));
    active_segments = all_segments(run_idx);
end

function target_values = resolve_target_values(target_T_values, t_sched, target_T_offsets)
    if isempty(target_T_values)
        target_values = t_sched + target_T_offsets;
    else
        target_values = target_T_values(:)';
    end
end

function [target_results, total_search_time_s, run_meta] = run_segment_target_search( ...
    active_segment, target_T_values, time_tolerance_s, control_refinement_steps, ...
    lambda_bounds, lambda_search_iters, lambda_seed_count, auto_adjust_dv, ...
    min_control_step_small_dx, max_lambda_expand_steps, min_lambda_value, max_lambda_value, ...
    project_root, Mass, gravity, Davis, inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

    fprintf('\n============================================================\n');
    fprintf('Loading data for %s (T_sched = %.0f s)...\n', ...
        active_segment.name, active_segment.T_sched);
    fprintf('Targets to search: [%s]\n', sprintf('%.0f ', target_T_values));

    route_file = fullfile(project_root, 'route', 'Guangxhou_line_7', active_segment.file);
    assert(exist(route_file, 'file') == 2, 'Route file tidak ditemukan: %s', route_file);
    data = load(route_file);

    if isfield(data, 'vel_profile')
        vp_raw = data.vel_profile;
    elseif isfield(data, 'speed_limit')
        vp_raw = data.speed_limit;
    else
        error('Data route %s tidak punya vel_profile/speed_limit.', active_segment.file);
    end
    if isfield(data, 'gradient')
        grad_raw = data.gradient;
    elseif isfield(data, 'slope')
        grad_raw = data.slope;
    else
        error('Data route %s tidak punya gradient/slope.', active_segment.file);
    end

    vel_lim = [vp_raw(:,1)*1000, vp_raw(:,2)/3.6];
    grad = [grad_raw(:,1)*1000, grad_raw(:,2)];
    total_dist = vel_lim(end,1);

    dx = 1;
    dv = 0.3;
    control_levels = [-1.0, -0.75, -0.5, -0.25, 0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0];
    a_max = 1.3;
    a_min = -1.3;
    epsilon_v = 0.5;

    dv_requested = dv;
    if auto_adjust_dv
        dv = recommend_velocity_grid(dx, dv_requested, a_max, max(vel_lim(:,2)));
        if abs(dv - dv_requested) > 1e-9
            fprintf('Auto-adjust dv from %.3f to %.3f for dx=%.1f to reduce speed-bin aliasing.\n', ...
                dv_requested, dv, dx);
        end
    end

    if dx <= 2
        original_steps = control_refinement_steps;
        control_refinement_steps = control_refinement_steps(control_refinement_steps >= min_control_step_small_dx - 1e-9);
        if numel(control_refinement_steps) < numel(original_steps)
            fprintf('Skipping control refinements finer than %.3f for dx=%.1f to keep runtime practical.\n', ...
                min_control_step_small_dx, dx);
        end
    end

    x_grid = 0 : dx : total_dist;
    v_max = max(vel_lim(:,2)) + 2;
    v_grid = 0 : dv : v_max;
    N_stages = length(x_grid);
    Nv = length(v_grid);
    Nc = length(control_levels);

    fprintf('Grid: %d stages, %d speeds, %d controls\n', N_stages, Nv, Nc);
    fprintf('Running adaptive target-time search...\n');

    n_targets = numel(target_T_values);
    overall_search_timer = tic;
    target_results = repmat(struct( ...
        'target_T', NaN, ...
        'T_total', NaN, ...
        'E_total_kWh', NaN, ...
        'time_error_s', NaN, ...
        'within_tolerance', false, ...
        'best_lambda', NaN, ...
        'control_step', NaN, ...
        'solve_time_s', NaN, ...
        'v_opt', [], ...
        't_opt', [], ...
        'x_opt', []), 1, n_targets);

    figure('Name', sprintf('DP %s target-time tracking', active_segment.name));
    hold on;
    colors = lines(max(n_targets, 1));
    plot_handles = gobjects(0);
    plot_labels = {};

    for idx = 1:n_targets
        target_T = target_T_values(idx);
        target_timer = tic;
        target_results(idx) = solve_for_target_time(target_T, time_tolerance_s, ...
            control_refinement_steps, lambda_bounds, lambda_search_iters, lambda_seed_count, ...
            max_lambda_expand_steps, min_lambda_value, max_lambda_value, ...
            vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
            Mass, gravity, Davis, inertial_mass_kg, ...
            Max_tractive_power, V1_traction, V2_traction, ...
            Max_brake_power, V1_brake, V2_brake);
        target_results(idx).solve_time_s = toc(target_timer);

        status_text = 'closest';
        if target_results(idx).within_tolerance
            status_text = 'matched';
        end

        fprintf(['Target = %6.1f s | Matched = %7.3f s | Error = %+7.3f s | ' ...
            'E = %6.3f kWh | lambda = %.5f | du = %.3f | solve = %.2f s | %s\n'], ...
            target_results(idx).target_T, target_results(idx).T_total, ...
            target_results(idx).time_error_s, target_results(idx).E_total_kWh, ...
            target_results(idx).best_lambda, target_results(idx).control_step, ...
            target_results(idx).solve_time_s, status_text);

        plot_handles(end+1) = plot(target_results(idx).x_opt/1000, ...
            target_results(idx).v_opt*3.6, 'LineWidth', 1.8, 'Color', colors(idx,:)); %#ok<SAGROW>
        plot_labels{end+1} = sprintf('T*=%.1fs -> %.3fs (E=%.3fkWh, du=%.3f)', ...
            target_results(idx).target_T, target_results(idx).T_total, ...
            target_results(idx).E_total_kWh, target_results(idx).control_step); %#ok<SAGROW>
    end

    limit_handle = stairs(vel_lim(:,1)/1000, vel_lim(:,2)*3.6, 'k--', 'LineWidth', 2);
    legend([plot_handles, limit_handle], [plot_labels, {'Speed Limit'}], 'Location', 'best');
    xlabel('Posisi (km)');
    ylabel('Kecepatan (km/h)');
    title(sprintf('Profil DP %s untuk target waktu yang diminta', active_segment.name));
    grid on;

    total_search_time_s = toc(overall_search_timer);
    run_meta = struct( ...
        'route_file', route_file, ...
        'dx', dx, ...
        'dv_requested', dv_requested, ...
        'dv_used', dv, ...
        'n_stages', N_stages, ...
        'n_speeds', Nv, ...
        'n_controls', Nc, ...
        'total_dist_m', total_dist);
end

function summary_text = build_segment_summary_text(active_segment, target_T_values, target_results, total_search_time_s, run_meta)
    lines = {
        sprintf('=== RINGKASAN TARGET TIME | %s ===', active_segment.name)
        sprintf('Route file            : %s', active_segment.file)
        sprintf('T_sched               : %.1f s', active_segment.T_sched)
        sprintf('Targets               : [%s]', sprintf('%.0f ', target_T_values))
        sprintf('Grid                  : dx=%.2f m | dv_req=%.2f m/s | dv_used=%.2f m/s | stages=%d | speeds=%d | controls=%d', ...
            run_meta.dx, run_meta.dv_requested, run_meta.dv_used, ...
            run_meta.n_stages, run_meta.n_speeds, run_meta.n_controls)
        sprintf('Distance              : %.1f m', run_meta.total_dist_m)
        ''
        sprintf('%-10s | %-10s | %-10s | %-12s | %-10s | %-10s | %-8s', ...
            'Target (s)', 'Match (s)', 'Error (s)', 'Energi (kWh)', 'Lambda', 'Solve (s)', 'Status')
        '-----------------------------------------------------------------------------------------'
        };

    for idx = 1:numel(target_results)
        status_text = 'closest';
        if target_results(idx).within_tolerance
            status_text = 'matched';
        end
        lines{end+1,1} = sprintf('%-10.1f | %-10.3f | %-+10.3f | %-12.3f | %-10.5f | %-10.2f | %-8s', ...
            target_results(idx).target_T, target_results(idx).T_total, ...
            target_results(idx).time_error_s, target_results(idx).E_total_kWh, ...
            target_results(idx).best_lambda, target_results(idx).solve_time_s, status_text); %#ok<AGROW>
    end

    matched_count = sum([target_results.within_tolerance]);
    avg_abs_error_s = mean(abs([target_results.time_error_s]), 'omitnan');
    avg_solve_time_s = mean([target_results.solve_time_s], 'omitnan');
    lines{end+1,1} = '-----------------------------------------------------------------------------------------';
    lines{end+1,1} = sprintf('Matched %d/%d | Avg |error| = %.3f s | Total search wall time = %.2f s | Avg solve = %.2f s', ...
        matched_count, numel(target_results), avg_abs_error_s, total_search_time_s, avg_solve_time_s);
    summary_text = sprintf('%s\n', lines{:});
end

function summary_text = build_overall_summary_text(segment_runs, overall_batch_time_s, summary_output_dir)
    lines = {
        '=== RINGKASAN SEMUA SEGMENT ==='
        sprintf('%-8s | %-10s | %-11s | %-8s | %-13s | %-12s', ...
            'Segment', 'T_sched', 'Targets', 'Match', 'Avg|Err| (s)', 'Solve (s)')
        '--------------------------------------------------------------------------------'
        };

    total_target_count = 0;
    total_matched_count = 0;
    for idx = 1:numel(segment_runs)
        target_results = segment_runs(idx).target_results;
        matched_count = sum([target_results.within_tolerance]);
        target_count = numel(target_results);
        avg_abs_error_s = mean(abs([target_results.time_error_s]), 'omitnan');
        target_label = format_target_list(segment_runs(idx).target_T_values);

        lines{end+1,1} = sprintf('%-8s | %-10.1f | %-11s | %-8s | %-13.3f | %-12.2f', ...
            segment_runs(idx).segment.name, segment_runs(idx).segment.T_sched, ...
            target_label, sprintf('%d/%d', matched_count, target_count), ...
            avg_abs_error_s, segment_runs(idx).total_search_time_s); %#ok<AGROW>

        total_target_count = total_target_count + target_count;
        total_matched_count = total_matched_count + matched_count;
    end

    lines{end+1,1} = '--------------------------------------------------------------------------------';
    lines{end+1,1} = sprintf('Matched total %d/%d | Overall wall time = %.2f s', ...
        total_matched_count, total_target_count, overall_batch_time_s);
    lines{end+1,1} = sprintf('Saved per-segment and aggregate summaries in: %s', summary_output_dir);
    summary_text = sprintf('%s\n', lines{:});
end

function target_label = format_target_list(target_values)
    if isempty(target_values)
        target_label = '-';
    elseif numel(target_values) == 1
        target_label = sprintf('%.0f', target_values);
    else
        target_label = sprintf('%.0f:%.0f', target_values(1), target_values(end));
    end
end

function write_text_file(file_path, file_text)
    fid = fopen(file_path, 'w');
    assert(fid ~= -1, 'Gagal membuka file summary untuk ditulis: %s', file_path);
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', file_text);
end

function result = solve_for_target_time(target_T, time_tolerance_s, ...
    control_refinement_steps, lambda_bounds, lambda_search_iters, lambda_seed_count, ...
    max_lambda_expand_steps, min_lambda_value, max_lambda_value, ...
    vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
    Mass, gravity, Davis, inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

    best_result = struct('target_T', target_T, 'T_total', NaN, 'E_total_kWh', NaN, ...
        'time_error_s', Inf, 'within_tolerance', false, 'best_lambda', NaN, ...
        'control_step', NaN, 'solve_time_s', NaN, 'v_opt', [], 't_opt', [], 'x_opt', []);

    for step = control_refinement_steps
        control_levels = build_control_levels(step);
        [lambda_low, lambda_high, low_candidate, high_candidate] = expand_lambda_bounds( ...
            target_T, step, control_levels, lambda_bounds, max_lambda_expand_steps, ...
            min_lambda_value, max_lambda_value, ...
            vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
            Mass, gravity, Davis, inertial_mass_kg, ...
            Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);

        if abs(lambda_high - lambda_low) <= 1e-12
            lambda_seed = lambda_low;
        elseif lambda_seed_count <= 2
            lambda_seed = unique([lambda_low, lambda_high]);
        else
            lambda_seed = unique(logspace(log10(lambda_low), log10(lambda_high), lambda_seed_count));
        end

        candidates = repmat(best_result, 1, numel(lambda_seed));

        for idx = 1:numel(lambda_seed)
            if abs(lambda_seed(idx) - lambda_low) <= 1e-12
                candidates(idx) = low_candidate;
            elseif abs(lambda_seed(idx) - lambda_high) <= 1e-12
                candidates(idx) = high_candidate;
            else
                candidates(idx) = evaluate_target_candidate(target_T, lambda_seed(idx), step, control_levels, ...
                    vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
                    Mass, gravity, Davis, inertial_mass_kg, ...
                    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
            end
            best_result = pick_better_result(best_result, candidates(idx), time_tolerance_s);
        end

        fprintf('  du=%.3f | lambda=[%.5f, %.5f] -> T=[%.2f, %.2f]\n', ...
            step, lambda_low, lambda_high, low_candidate.T_total, high_candidate.T_total);

        if best_result.within_tolerance
            result = best_result;
            return;
        end

        [has_bracket, lam_lo, lam_hi, res_lo, res_hi] = find_lambda_bracket(candidates, target_T);
        if ~has_bracket
            continue;
        end

        for iter = 1:lambda_search_iters
            lam_mid = 0.5 * (lam_lo + lam_hi);
            res_mid = evaluate_target_candidate(target_T, lam_mid, step, control_levels, ...
                vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
                Mass, gravity, Davis, inertial_mass_kg, ...
                Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
            best_result = pick_better_result(best_result, res_mid, time_tolerance_s);

            if best_result.within_tolerance
                result = best_result;
                return;
            end

            if abs(res_mid.time_error_s) <= 1e-9 || abs(lam_hi - lam_lo) <= 1e-6
                break;
            end

            if res_mid.T_total > target_T
                lam_lo = lam_mid;
                res_lo = res_mid;
            else
                lam_hi = lam_mid;
                res_hi = res_mid;
            end

            if abs(res_lo.T_total - res_hi.T_total) <= 1e-6
                break;
            end
        end
    end

    result = best_result;
end

function candidate = evaluate_target_candidate(target_T, lambda_value, control_step, control_levels, ...
    vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
    Mass, gravity, Davis, inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

    [T_total, E_total_kWh, v_opt, t_opt, x_opt] = solve_dp(lambda_value, ...
        vel_lim, grad, x_grid, v_grid, control_levels, dx, a_max, a_min, epsilon_v, ...
        Mass, gravity, Davis, inertial_mass_kg, ...
        Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);

    candidate = struct( ...
        'target_T', target_T, ...
        'T_total', T_total, ...
        'E_total_kWh', E_total_kWh, ...
        'time_error_s', T_total - target_T, ...
        'within_tolerance', false, ...
        'best_lambda', lambda_value, ...
        'control_step', control_step, ...
        'solve_time_s', NaN, ...
        'v_opt', v_opt, ...
        't_opt', t_opt, ...
        'x_opt', x_opt);
end

function best_result = pick_better_result(best_result, candidate, time_tolerance_s)
    candidate.within_tolerance = abs(candidate.time_error_s) <= time_tolerance_s;
    best_result.within_tolerance = abs(best_result.time_error_s) <= time_tolerance_s;

    if candidate.within_tolerance && ~best_result.within_tolerance
        best_result = candidate;
        return;
    end

    if candidate.within_tolerance && best_result.within_tolerance
        if candidate.E_total_kWh < best_result.E_total_kWh
            best_result = candidate;
        end
        return;
    end

    if abs(candidate.time_error_s) < abs(best_result.time_error_s) - 1e-9
        best_result = candidate;
    elseif abs(abs(candidate.time_error_s) - abs(best_result.time_error_s)) <= 1e-9 && ...
            candidate.E_total_kWh < best_result.E_total_kWh
        best_result = candidate;
    end
end

function [has_bracket, lam_lo, lam_hi, res_lo, res_hi] = find_lambda_bracket(candidates, target_T)
    has_bracket = false;
    lam_lo = NaN;
    lam_hi = NaN;
    res_lo = struct();
    res_hi = struct();

    [~, order] = sort([candidates.best_lambda]);
    candidates = candidates(order);

    for idx = 1 : numel(candidates) - 1
        T_left = candidates(idx).T_total;
        T_right = candidates(idx + 1).T_total;
        if T_left >= target_T && T_right <= target_T
            has_bracket = true;
            lam_lo = candidates(idx).best_lambda;
            lam_hi = candidates(idx + 1).best_lambda;
            res_lo = candidates(idx);
            res_hi = candidates(idx + 1);
            return;
        end
    end
end

function iv_match = nearest_feasible_policy_bin(policy_column, v_current, v_grid)
    feasible_bins = find(isfinite(policy_column));
    if isempty(feasible_bins)
        iv_match = 0;
        return;
    end

    [~, local_idx] = min(abs(v_grid(feasible_bins) - v_current));
    iv_match = feasible_bins(local_idx);
end

function control_levels = build_control_levels(step)
    raw_levels = -1 : step : 1;
    control_levels = unique(round(raw_levels / step) * step);
    if control_levels(1) > -1
        control_levels = [-1, control_levels];
    end
    if control_levels(end) < 1
        control_levels = [control_levels, 1];
    end
end

function dv_used = recommend_velocity_grid(dx, dv_requested, a_max, v_lim_max)
    v_ref = max(v_lim_max / 2, 5);
    dv_raw = a_max * dx / v_ref;
    dv_target = max(0.05, round(dv_raw / 0.05) * 0.05);
    dv_used = min(dv_requested, dv_target);
end

function J_next = interpolate_cost_to_go(v_next, v_grid, J_column)
    if v_next <= v_grid(1)
        J_next = J_column(1);
        return;
    end

    if v_next >= v_grid(end)
        J_next = J_column(end);
        return;
    end

    iv_hi = find(v_grid >= v_next, 1, 'first');
    iv_lo = iv_hi - 1;
    J_lo = J_column(iv_lo);
    J_hi = J_column(iv_hi);

    if isinf(J_lo) && isinf(J_hi)
        J_next = inf;
        return;
    end

    if isinf(J_lo)
        J_next = J_hi;
        return;
    end

    if isinf(J_hi)
        J_next = J_lo;
        return;
    end

    alpha = (v_next - v_grid(iv_lo)) / (v_grid(iv_hi) - v_grid(iv_lo));
    J_next = (1 - alpha) * J_lo + alpha * J_hi;
end

function [lambda_low, lambda_high, low_candidate, high_candidate] = expand_lambda_bounds( ...
    target_T, control_step, control_levels, lambda_bounds, max_lambda_expand_steps, ...
    min_lambda_value, max_lambda_value, ...
    vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
    Mass, gravity, Davis, inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

    lambda_low = max(lambda_bounds(1), min_lambda_value);
    lambda_high = min(lambda_bounds(2), max_lambda_value);

    low_candidate = evaluate_target_candidate(target_T, lambda_low, control_step, control_levels, ...
        vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
        Mass, gravity, Davis, inertial_mass_kg, ...
        Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
    high_candidate = evaluate_target_candidate(target_T, lambda_high, control_step, control_levels, ...
        vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
        Mass, gravity, Davis, inertial_mass_kg, ...
        Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);

    expand_iter = 0;
    while low_candidate.T_total < target_T && lambda_low > min_lambda_value + 1e-12 && expand_iter < max_lambda_expand_steps
        lambda_low = max(min_lambda_value, lambda_low / 2);
        low_candidate = evaluate_target_candidate(target_T, lambda_low, control_step, control_levels, ...
            vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
            Mass, gravity, Davis, inertial_mass_kg, ...
            Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
        expand_iter = expand_iter + 1;
    end

    expand_iter = 0;
    while high_candidate.T_total > target_T && lambda_high < max_lambda_value - 1e-12 && expand_iter < max_lambda_expand_steps
        lambda_high = min(max_lambda_value, lambda_high * 2);
        high_candidate = evaluate_target_candidate(target_T, lambda_high, control_step, control_levels, ...
            vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
            Mass, gravity, Davis, inertial_mass_kg, ...
            Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
        expand_iter = expand_iter + 1;
    end
end

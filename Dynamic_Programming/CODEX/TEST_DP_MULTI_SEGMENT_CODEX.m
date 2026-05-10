%% =========================================================================
% TEST_DP_MULTI_SEGMENT_CODEX.m
% Improved DP target-time search:
% - keeps searching after the first time-feasible match to reduce energy
% - applies local lambda polishing around matched solutions
% - uses denser near-coast control levels for finer energy shaping
% - saves smoother plots and summaries under Dynamic_Programming/CODEX/results
% =========================================================================
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(genpath(project_root));
route_direction = 'up';

ALL_SEGMENTS = get_guangzhou_line7_catalog(route_direction);

selected_segments = 'IS01';  % bisa juga 'all', 1..8, [1 3 5], {'IS01','IS03'}
segments_to_run = resolve_selected_segments(selected_segments, ALL_SEGMENTS);

target_T_offsets = 0;    % quick sanity run; set 0:4 for [T_sched, T_sched+1, ...]
target_T_values = [];    % isi manual kalau mau override, contoh: [129, 130, 131]
time_tolerance_s = 0.5;
control_refinement_steps = [0.15];
lambda_bounds = [0.005, 0.2];
lambda_search_iters = 12;
lambda_seed_count = 2;
auto_adjust_dv = true;
min_control_step_small_dx = 0.05;
max_lambda_expand_steps = 10;
min_lambda_value = 1e-4;
max_lambda_value = 10;

improve = struct();
improve.continue_search_after_first_match = true;
improve.refine_lambda_relative_span = 0.20;
improve.enable_energy_polish = true;
improve.energy_polish_relative_span = 0.08;
improve.energy_polish_lambda_samples = 5;

display_cfg = struct();
display_cfg.show_plots = false;
display_cfg.save_plots = true;
display_cfg.show_raw_points = false;
display_cfg.smoothing_factor = 8;

results_output_dir = fullfile(script_dir, 'results');
if exist(results_output_dir, 'dir') ~= 7
	mkdir(results_output_dir);
end

fprintf('Loading rolling stock...\n');
run(fullfile(project_root, 'rollingstocks', 'rollingstock_Guangzhou_L7.m'));
inertial_mass_kg = Mass * (1 + lambda) * 1000;

overall_batch_timer = tic;
segment_runs = repmat(struct( ...
	'segment', struct('name', '', 'file', '', 'T_sched', NaN, 'direction', ''), ...
	'target_T_values', [], ...
	'target_results', [], ...
	'total_search_time_s', NaN, ...
	'summary_file', '', ...
	'results_file', '', ...
	'plot_files', struct('speed_profile', '', 'energy_curve', ''), ...
	'run_meta', struct()), 1, numel(segments_to_run));

for seg_run_idx = 1:numel(segments_to_run)
	active_segment = segments_to_run(seg_run_idx);
	segment_target_T_values = resolve_target_values(target_T_values, active_segment.T_sched, target_T_offsets);

	[target_results, total_search_time_s, run_meta] = run_segment_target_search( ...
		active_segment, segment_target_T_values, time_tolerance_s, ...
		control_refinement_steps, lambda_bounds, lambda_search_iters, ...
		lambda_seed_count, auto_adjust_dv, min_control_step_small_dx, ...
		max_lambda_expand_steps, min_lambda_value, max_lambda_value, ...
		improve, display_cfg, project_root, results_output_dir, ...
		Mass, gravity, Davis, inertial_mass_kg, ...
		Max_tractive_power, V1_traction, V2_traction, ...
		Max_brake_power, V1_brake, V2_brake);

	segment_summary_text = build_segment_summary_text( ...
		active_segment, segment_target_T_values, target_results, total_search_time_s, run_meta);
	fprintf('\n%s', segment_summary_text);

	route_base_name = erase(active_segment.file, '.mat');
	segment_summary_file = fullfile(results_output_dir, sprintf('%s_codex_summary.txt', route_base_name));
	segment_results_file = fullfile(results_output_dir, sprintf('%s_codex_results.mat', route_base_name));
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
		'plot_files', run_meta.plot_files, ...
		'run_meta', run_meta);
end

overall_batch_time_s = toc(overall_batch_timer);
aggregate_summary_text = build_overall_summary_text(segment_runs, overall_batch_time_s, results_output_dir);
fprintf('\n%s', aggregate_summary_text);

aggregate_summary_file = fullfile(results_output_dir, 'all_segments_codex_summary.txt');
aggregate_results_file = fullfile(results_output_dir, 'all_segments_codex_results.mat');
write_text_file(aggregate_summary_file, aggregate_summary_text);
save(aggregate_results_file, 'segment_runs', 'overall_batch_time_s', 'segments_to_run', ...
	'target_T_offsets', 'target_T_values', 'improve', 'display_cfg');

fprintf('CODEX results saved in: %s\n', results_output_dir);

function F = traction_force(v, Pmax, V1, V2)
	if v <= V1
		F = Pmax / V1;
	elseif v <= V2
		F = Pmax / v;
	else
		F = Pmax * V2 / v^2;
	end
end

function F = brake_force(v, Pmax, V1, V2)
	if v <= V1
		F = Pmax / V1;
	elseif v <= V2
		F = Pmax / v;
	else
		F = Pmax * V2 / v^2;
	end
end

function [T_total, E_total_kWh, v_opt, t_opt, x_opt] = solve_dp(lambda_value, ...
	vel_lim, grad, x_grid, v_grid, control_levels, dx, a_max, a_min, epsilon_v, ...
	Mass, gravity, Davis, inertial_mass_kg, ...
	Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

	N_stages = length(x_grid);
	Nv = length(v_grid);
	Nc = length(control_levels);

	J = inf(Nv, N_stages);
	policy = nan(Nv, N_stages - 1);

	for iv = 1:Nv
		if v_grid(iv) <= epsilon_v
			J(iv, end) = 0;
		end
	end

	for k = N_stages - 1 : -1 : 1
		x_cur = x_grid(k);
		v_lim = interp1(vel_lim(:,1), vel_lim(:,2), x_cur, 'previous', vel_lim(end,2));
		v_lim_next = interp1(vel_lim(:,1), vel_lim(:,2), x_grid(k + 1), 'previous', vel_lim(end,2));
		grad_permil = interp1(grad(:,1), grad(:,2), x_cur, 'previous', 0);
		theta = atan(grad_permil / 1000);
		F_grade = Mass * 1000 * gravity * sin(theta);

		for iv = 1:Nv
			v = v_grid(iv);
			if v > v_lim + 1e-6
				continue;
			end

			R_davis = (Davis(1) + Davis(2) * v + Davis(3) * v^2) * 1000;

			for ic = 1:Nc
				c = control_levels(ic);
				if c >= 0
					F_max = traction_force(v, Max_tractive_power, V1_traction, V2_traction);
					F_applied = c * F_max;
				else
					F_max = brake_force(v, Max_brake_power, V1_brake, V2_brake);
					F_applied = c * F_max;
				end

				F_net = F_applied - R_davis - F_grade;
				a = F_net / inertial_mass_kg;
				if a > a_max || a < a_min
					continue;
				end

				v_sq = v^2 + 2 * a * dx;
				if v_sq < 0
					continue;
				end
				v_next = sqrt(v_sq);
				if v_next > v_lim_next + 1e-6
					continue;
				end

				dt = 2 * dx / (v + max(v_next, 0.01));
				dE_kWh = max(0, F_applied) * dx / 3.6e6;

				if k == N_stages - 1
					if v_next <= epsilon_v
						J_cand = dE_kWh + lambda_value * dt;
						if J_cand < J(iv, k)
							J(iv, k) = J_cand;
							policy(iv, k) = c;
						end
					end
				else
					J_next = interpolate_cost_to_go(v_next, v_grid, J(:, k + 1));
					if isinf(J_next)
						continue;
					end
					J_cand = dE_kWh + lambda_value * dt + J_next;
					if J_cand < J(iv, k)
						J(iv, k) = J_cand;
						policy(iv, k) = c;
					end
				end
			end

			if k == N_stages - 1 && v > epsilon_v
				a_stop = -v^2 / (2 * dx);
				if a_stop >= a_min - 1e-9
					F_brake_cap = brake_force(v, Max_brake_power, V1_brake, V2_brake);
					F_applied_stop = inertial_mass_kg * a_stop + R_davis + F_grade;
					if F_applied_stop <= 0 && F_applied_stop >= -F_brake_cap - 1e-9
						c_stop = max(-1, F_applied_stop / max(F_brake_cap, 1));
						dt_stop = 2 * dx / max(v, 0.01);
						J_cand = lambda_value * dt_stop;
						if J_cand < J(iv, k)
							J(iv, k) = J_cand;
							policy(iv, k) = c_stop;
						end
					end
				end
			end
		end
	end

	[~, iv0] = min(abs(v_grid));
	if isinf(J(iv0, 1))
		warning('solve_dp: J(v=0, stage=1) = inf. No feasible path found (lambda=%.4g).', lambda_value);
	end

	iv = iv0;
	v_opt = zeros(1, N_stages);
	t_opt = zeros(1, N_stages);
	E_opt = zeros(1, N_stages);
	x_opt = x_grid;

	for k = 1:N_stages - 1
		c = policy(iv, k);
		if isnan(c)
			iv_fallback = nearest_feasible_policy_bin(policy(:, k), v_opt(k), v_grid);
			if iv_fallback == 0
				warning('solve_dp: no feasible control at stage %d (iv=%d). Stopping.', k, iv);
				break;
			end
			c = policy(iv_fallback, k);
		end

		v = v_opt(k);
		x_cur = x_grid(k);
		grad_permil = interp1(grad(:,1), grad(:,2), x_cur, 'previous', 0);
		theta = atan(grad_permil / 1000);
		F_grade = Mass * 1000 * gravity * sin(theta);
		R_davis = (Davis(1) + Davis(2) * v + Davis(3) * v^2) * 1000;

		if c >= 0
			F_max = traction_force(v, Max_tractive_power, V1_traction, V2_traction);
			F_applied = c * F_max;
		else
			F_max = brake_force(v, Max_brake_power, V1_brake, V2_brake);
			F_applied = c * F_max;
		end

		a = (F_applied - R_davis - F_grade) / inertial_mass_kg;
		v_sq = v^2 + 2 * a * dx;
		v_next = sqrt(max(v_sq, 0));

		if k == N_stages - 1
			v_next = min(v_next, epsilon_v);
		end

		[~, iv] = min(abs(v_grid - v_next));
		v_opt(k + 1) = v_next;
		dt = 2 * dx / (v + max(v_next, 0.01));
		t_opt(k + 1) = t_opt(k) + dt;
		E_opt(k + 1) = E_opt(k) + max(0, F_applied) * dx / 3.6e6;
	end

	T_total = t_opt(end);
	E_total_kWh = E_opt(end);
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
	improve, display_cfg, project_root, results_output_dir, ...
	Mass, gravity, Davis, inertial_mass_kg, ...
	Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

	fprintf('\n============================================================\n');
	fprintf('Loading data for %s (T_sched = %.0f s)...\n', ...
		active_segment.name, active_segment.T_sched);
	fprintf('Targets to search: [%s]\n', sprintf('%.0f ', target_T_values));

	route_file = resolve_project_route_path(project_root, active_segment.file, active_segment.direction);
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

	vel_lim = [vp_raw(:,1) * 1000, vp_raw(:,2) / 3.6];
	grad = [grad_raw(:,1) * 1000, grad_raw(:,2)];
	total_dist = vel_lim(end,1);

	dx = 1;
	dv = 0.3;
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
	base_control_levels = build_control_levels(control_refinement_steps(1));
	N_stages = length(x_grid);
	Nv = length(v_grid);
	Nc = numel(base_control_levels);

	fprintf('Grid: %d stages, %d speeds, %d base-controls\n', N_stages, Nv, Nc);
	fprintf('Running energy-polished target-time search...\n');

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

	for idx = 1:n_targets
		target_T = target_T_values(idx);
		target_timer = tic;
		target_results(idx) = solve_for_target_time( ...
			target_T, time_tolerance_s, control_refinement_steps, ...
			lambda_bounds, lambda_search_iters, lambda_seed_count, ...
			max_lambda_expand_steps, min_lambda_value, max_lambda_value, ...
			improve, vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
			Mass, gravity, Davis, inertial_mass_kg, ...
			Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
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
	end

	total_search_time_s = toc(overall_search_timer);
	plot_files = save_segment_plots(active_segment, vel_lim, target_results, display_cfg, results_output_dir);

	run_meta = struct( ...
		'route_file', route_file, ...
		'dx', dx, ...
		'dv_requested', dv_requested, ...
		'dv_used', dv, ...
		'n_stages', N_stages, ...
		'n_speeds', Nv, ...
		'n_controls', Nc, ...
		'total_dist_m', total_dist, ...
		'plot_files', plot_files, ...
		'improve', improve);
end

function summary_text = build_segment_summary_text(active_segment, target_T_values, target_results, total_search_time_s, run_meta)
	lines = {
		sprintf('=== CODEX TARGET-TIME SUMMARY | %s ===', active_segment.name)
		sprintf('Route file            : %s', active_segment.file)
		sprintf('T_sched               : %.1f s', active_segment.T_sched)
		sprintf('Targets               : [%s]', sprintf('%.0f ', target_T_values))
		sprintf('Grid                  : dx=%.2f m | dv_req=%.2f m/s | dv_used=%.2f m/s | stages=%d | speeds=%d | controls=%d', ...
			run_meta.dx, run_meta.dv_requested, run_meta.dv_used, ...
			run_meta.n_stages, run_meta.n_speeds, run_meta.n_controls)
		sprintf('Distance              : %.1f m', run_meta.total_dist_m)
		sprintf('Energy polish         : %d | continue-after-match: %d', ...
			run_meta.improve.enable_energy_polish, run_meta.improve.continue_search_after_first_match)
		sprintf('Speed plot            : %s', run_meta.plot_files.speed_profile)
		sprintf('Energy plot           : %s', run_meta.plot_files.energy_curve)
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
	best_energy_kWh = min([target_results.E_total_kWh]);

	lines{end+1,1} = '-----------------------------------------------------------------------------------------';
	lines{end+1,1} = sprintf('Matched %d/%d | Avg |error| = %.3f s | Best energy = %.3f kWh | Total search wall time = %.2f s | Avg solve = %.2f s', ...
		matched_count, numel(target_results), avg_abs_error_s, best_energy_kWh, total_search_time_s, avg_solve_time_s);
	summary_text = sprintf('%s\n', lines{:});
end

function summary_text = build_overall_summary_text(segment_runs, overall_batch_time_s, results_output_dir)
	lines = {
		'=== CODEX SUMMARY FOR ALL SEGMENTS ==='
		sprintf('%-8s | %-10s | %-11s | %-8s | %-13s | %-12s | %-10s', ...
			'Segment', 'T_sched', 'Targets', 'Match', 'Avg|Err| (s)', 'Solve (s)', 'Best E')
		'-------------------------------------------------------------------------------------------'
		};

	total_target_count = 0;
	total_matched_count = 0;
	for idx = 1:numel(segment_runs)
		target_results = segment_runs(idx).target_results;
		matched_count = sum([target_results.within_tolerance]);
		target_count = numel(target_results);
		avg_abs_error_s = mean(abs([target_results.time_error_s]), 'omitnan');
		target_label = format_target_list(segment_runs(idx).target_T_values);
		best_energy_kWh = min([target_results.E_total_kWh]);

		lines{end+1,1} = sprintf('%-8s | %-10.1f | %-11s | %-8s | %-13.3f | %-12.2f | %-10.3f', ...
			segment_runs(idx).segment.name, segment_runs(idx).segment.T_sched, ...
			target_label, sprintf('%d/%d', matched_count, target_count), ...
			avg_abs_error_s, segment_runs(idx).total_search_time_s, best_energy_kWh); %#ok<AGROW>

		total_target_count = total_target_count + target_count;
		total_matched_count = total_matched_count + matched_count;
	end

	lines{end+1,1} = '-------------------------------------------------------------------------------------------';
	lines{end+1,1} = sprintf('Matched total %d/%d | Overall wall time = %.2f s', ...
		total_matched_count, total_target_count, overall_batch_time_s);
	lines{end+1,1} = sprintf('Saved summaries, MAT files, and plots in: %s', results_output_dir);
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
	cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
	fprintf(fid, '%s', file_text);
end

function plot_files = save_segment_plots(active_segment, vel_lim, target_results, display_cfg, results_output_dir)
	route_base_name = erase(active_segment.file, '.mat');
	speed_plot_file = fullfile(results_output_dir, sprintf('%s_codex_speed_profile.png', route_base_name));
	energy_plot_file = fullfile(results_output_dir, sprintf('%s_codex_energy_curve.png', route_base_name));

	figure_visibility = 'off';
	if display_cfg.show_plots
		figure_visibility = 'on';
	end

	speed_fig = figure('Name', sprintf('CODEX %s speed profile', active_segment.name), ...
		'Visible', figure_visibility, 'Color', 'w');
	hold on;
	colors = lines(max(numel(target_results), 1));
	plot_handles = gobjects(0);
	plot_labels = {};

	for idx = 1:numel(target_results)
		[x_smooth, v_smooth] = smooth_profile_curve( ...
			target_results(idx).x_opt, target_results(idx).v_opt * 3.6, ...
			display_cfg.smoothing_factor, vel_lim);
		plot_handles(end + 1) = plot(x_smooth / 1000, v_smooth, ...
			'LineWidth', 2.0, 'Color', colors(idx,:)); %#ok<SAGROW>
		plot_labels{end + 1} = sprintf('T*=%.1fs -> %.3fs (E=%.3fkWh)', ...
			target_results(idx).target_T, target_results(idx).T_total, ...
			target_results(idx).E_total_kWh); %#ok<SAGROW>

		if display_cfg.show_raw_points
			plot(target_results(idx).x_opt / 1000, target_results(idx).v_opt * 3.6, ...
				':', 'Color', colors(idx,:), 'LineWidth', 0.8, 'HandleVisibility', 'off');
		end
	end

	limit_handle = stairs(vel_lim(:,1) / 1000, vel_lim(:,2) * 3.6, 'k--', 'LineWidth', 2);
	legend([plot_handles, limit_handle], [plot_labels, {'Speed Limit'}], 'Location', 'best');
	xlabel('Posisi (km)');
	ylabel('Kecepatan (km/h)');
	title(sprintf('CODEX Smoothed DP Speed Profile | %s', active_segment.name));
	grid on;
	exportgraphics(speed_fig, speed_plot_file, 'Resolution', 180);
	if ~display_cfg.show_plots
		close(speed_fig);
	end

	energy_fig = figure('Name', sprintf('CODEX %s energy curve', active_segment.name), ...
		'Visible', figure_visibility, 'Color', 'w');
	hold on;
	target_vector = [target_results.target_T];
	energy_vector = [target_results.E_total_kWh];
	plot(target_vector, energy_vector, '-o', 'LineWidth', 2.0, ...
		'MarkerFaceColor', [0.1 0.45 0.75], 'Color', [0.1 0.45 0.75]);
	xlabel('Target time (s)');
	ylabel('Traction energy (kWh)');
	title(sprintf('CODEX Energy vs Target Time | %s', active_segment.name));
	grid on;
	exportgraphics(energy_fig, energy_plot_file, 'Resolution', 180);
	if ~display_cfg.show_plots
		close(energy_fig);
	end

	plot_files = struct('speed_profile', speed_plot_file, 'energy_curve', energy_plot_file);
end

function [x_smooth, v_smooth] = smooth_profile_curve(x_raw, v_raw_kmh, smoothing_factor, vel_lim)
	x_raw = x_raw(:)';
	v_raw_kmh = v_raw_kmh(:)';
	valid_mask = [true, diff(x_raw) > 0];
	x_valid = x_raw(valid_mask);
	v_valid_kmh = v_raw_kmh(valid_mask);

	if numel(x_valid) < 3
		x_smooth = x_valid;
		v_smooth = max(v_valid_kmh, 0);
		return;
	end

	n_query = max(numel(x_valid), smoothing_factor * numel(x_valid));
	x_smooth = linspace(x_valid(1), x_valid(end), n_query);
	v_smooth = pchip(x_valid, v_valid_kmh, x_smooth);
	v_lim_kmh = interp1(vel_lim(:,1), vel_lim(:,2) * 3.6, x_smooth, 'previous', vel_lim(end,2) * 3.6);
	v_smooth = min(max(v_smooth, 0), v_lim_kmh);
end

function result = solve_for_target_time(target_T, time_tolerance_s, ...
	control_refinement_steps, lambda_bounds, lambda_search_iters, lambda_seed_count, ...
	max_lambda_expand_steps, min_lambda_value, max_lambda_value, improve, ...
	vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
	Mass, gravity, Davis, inertial_mass_kg, ...
	Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

	best_result = struct('target_T', target_T, 'T_total', NaN, 'E_total_kWh', NaN, ...
		'time_error_s', Inf, 'within_tolerance', false, 'best_lambda', NaN, ...
		'control_step', NaN, 'solve_time_s', NaN, 'v_opt', [], 't_opt', [], 'x_opt', []);

	for step = control_refinement_steps
		control_levels = build_control_levels(step);
		step_lambda_bounds = lambda_bounds;
		if isfinite(best_result.best_lambda)
			step_lambda_bounds = [ ...
				max(min_lambda_value, best_result.best_lambda * (1 - improve.refine_lambda_relative_span)), ...
				min(max_lambda_value, best_result.best_lambda * (1 + improve.refine_lambda_relative_span))];
		end

		[lambda_low, lambda_high, low_candidate, high_candidate] = expand_lambda_bounds( ...
			target_T, step, control_levels, step_lambda_bounds, max_lambda_expand_steps, ...
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
				candidates(idx) = evaluate_target_candidate( ...
					target_T, lambda_seed(idx), step, control_levels, ...
					vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
					Mass, gravity, Davis, inertial_mass_kg, ...
					Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
			end
			best_result = pick_better_result(best_result, candidates(idx), time_tolerance_s);
		end

		fprintf('  du=%.3f | lambda=[%.5f, %.5f] -> T=[%.2f, %.2f]\n', ...
			step, lambda_low, lambda_high, low_candidate.T_total, high_candidate.T_total);

		[has_bracket, lam_lo, lam_hi, res_lo, res_hi] = find_lambda_bracket(candidates, target_T);
		if has_bracket
			for iter = 1:lambda_search_iters
				lam_mid = 0.5 * (lam_lo + lam_hi);
				res_mid = evaluate_target_candidate( ...
					target_T, lam_mid, step, control_levels, ...
					vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
					Mass, gravity, Davis, inertial_mass_kg, ...
					Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
				best_result = pick_better_result(best_result, res_mid, time_tolerance_s);

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

		if best_result.within_tolerance && ~improve.continue_search_after_first_match
			break;
		end
	end

	if improve.enable_energy_polish && best_result.within_tolerance
		best_result = polish_matched_solution( ...
			best_result, target_T, time_tolerance_s, control_refinement_steps, ...
			improve.energy_polish_relative_span, improve.energy_polish_lambda_samples, ...
			min_lambda_value, max_lambda_value, ...
			vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
			Mass, gravity, Davis, inertial_mass_kg, ...
			Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
	end

	result = best_result;
end

function best_result = polish_matched_solution(best_result, target_T, time_tolerance_s, ...
	control_refinement_steps, energy_polish_relative_span, energy_polish_lambda_samples, ...
	min_lambda_value, max_lambda_value, ...
	vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
	Mass, gravity, Davis, inertial_mass_kg, ...
	Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

	candidate_steps = unique(control_refinement_steps(control_refinement_steps <= best_result.control_step + 1e-9));
	if isempty(candidate_steps)
		candidate_steps = best_result.control_step;
	end

	lambda_low = max(min_lambda_value, best_result.best_lambda * (1 - energy_polish_relative_span));
	lambda_high = min(max_lambda_value, best_result.best_lambda * (1 + energy_polish_relative_span));
	lambda_values = unique([ ...
		linspace(lambda_low, lambda_high, energy_polish_lambda_samples), ...
		best_result.best_lambda]);

	for step = candidate_steps
		control_levels = build_control_levels(step);
		for idx = 1:numel(lambda_values)
			candidate = evaluate_target_candidate( ...
				target_T, lambda_values(idx), step, control_levels, ...
				vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
				Mass, gravity, Davis, inertial_mass_kg, ...
				Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
			best_result = pick_better_result(best_result, candidate, time_tolerance_s);
		end
	end
end

function candidate = evaluate_target_candidate(target_T, lambda_value, control_step, control_levels, ...
	vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
	Mass, gravity, Davis, inertial_mass_kg, ...
	Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake)

	[T_total, E_total_kWh, v_opt, t_opt, x_opt] = solve_dp( ...
		lambda_value, vel_lim, grad, x_grid, v_grid, control_levels, ...
		dx, a_max, a_min, epsilon_v, Mass, gravity, Davis, inertial_mass_kg, ...
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
		if candidate.E_total_kWh < best_result.E_total_kWh - 1e-9
			best_result = candidate;
		elseif abs(candidate.E_total_kWh - best_result.E_total_kWh) <= 1e-9 && ...
				abs(candidate.time_error_s) < abs(best_result.time_error_s)
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

	for idx = 1:numel(candidates) - 1
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
	control_levels = -1 : step : 1;
	if step <= 0.05
		control_levels = [control_levels, -0.12, -0.08, -0.04, -0.02, 0, 0.02, 0.04, 0.08, 0.12];
	else
		control_levels = [control_levels, -0.10, -0.05, 0, 0.05, 0.10];
	end

	control_levels = unique(max(-1, min(1, control_levels)));
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

	low_candidate = evaluate_target_candidate( ...
		target_T, lambda_low, control_step, control_levels, ...
		vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
		Mass, gravity, Davis, inertial_mass_kg, ...
		Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
	high_candidate = evaluate_target_candidate( ...
		target_T, lambda_high, control_step, control_levels, ...
		vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
		Mass, gravity, Davis, inertial_mass_kg, ...
		Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);

	expand_iter = 0;
	while low_candidate.T_total < target_T && lambda_low > min_lambda_value + 1e-12 && expand_iter < max_lambda_expand_steps
		lambda_low = max(min_lambda_value, lambda_low / 2);
		low_candidate = evaluate_target_candidate( ...
			target_T, lambda_low, control_step, control_levels, ...
			vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
			Mass, gravity, Davis, inertial_mass_kg, ...
			Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
		expand_iter = expand_iter + 1;
	end

	expand_iter = 0;
	while high_candidate.T_total > target_T && lambda_high < max_lambda_value - 1e-12 && expand_iter < max_lambda_expand_steps
		lambda_high = min(max_lambda_value, lambda_high * 2);
		high_candidate = evaluate_target_candidate( ...
			target_T, lambda_high, control_step, control_levels, ...
			vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
			Mass, gravity, Davis, inertial_mass_kg, ...
			Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake);
		expand_iter = expand_iter + 1;
	end
end

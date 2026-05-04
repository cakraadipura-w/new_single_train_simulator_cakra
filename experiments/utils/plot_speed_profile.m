function plot_speed_profile(mat_file, ri, pi_idx)
% PLOT_SPEED_PROFILE  Plot speed vs distance for a stored Pareto solution.
%
% Usage:
%   plot_speed_profile('experiment_results/IS04/E1_A_vs_D.mat')
%       Auto-selects: improved-strategy run with largest Pareto front,
%       midpoint of that front.
%
%   plot_speed_profile('experiment_results/IS04/E1_A_vs_D.mat', 15)
%       Use results(15), auto-select midpoint.
%
%   plot_speed_profile('experiment_results/IS04/E1_A_vs_D.mat', 15, 3)
%       Use results(15), Pareto point index 3 (sorted ascending by T).
%
% The .mat file must be produced by run_E1_ablation.m and contain a
% 'results' struct with fields: X, F, strategy, use_improved,
% route_file, rs_file, config_id, config_desc.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

if nargin < 1 || isempty(mat_file)
    error('Provide path to a .mat file (e.g. E1_A_vs_D.mat)');
end

data = load(mat_file);
if ~isfield(data, 'results')
    error('No ''results'' struct found in %s', mat_file);
end
results = data.results;
n_res   = numel(results);

%% --- Print available configs ---
fprintf('\n=== Available solutions in: %s ===\n', mat_file);
cfg_ids = unique({results.config_id}, 'stable');
for ci = 1:numel(cfg_ids)
    mask     = strcmp({results.config_id}, cfg_ids{ci});
    ri_list  = find(mask);
    pf_sizes = arrayfun(@(r) size(r.F,1), results(mask));
    desc     = results(ri_list(1)).config_desc;
    strat    = results(ri_list(1)).strategy;
    fprintf('  Config %-3s [%s, %s]: %d runs | PF size mean=%.1f max=%d\n', ...
        cfg_ids{ci}, desc, strat, numel(ri_list), mean(pf_sizes), max(pf_sizes));
end
fprintf('\n');

%% --- Auto-select result index ri ---
if nargin < 2 || isempty(ri)
    % Prefer 'improved' strategy; among those pick run with largest Pareto front
    imp_mask = false(1, n_res);
    for k = 1:n_res
        if isfield(results(k),'use_improved') && results(k).use_improved
            imp_mask(k) = true;
        elseif isfield(results(k),'strategy') && strcmp(results(k).strategy,'improved')
            imp_mask(k) = true;
        end
    end
    if ~any(imp_mask)
        imp_mask = true(1, n_res);  % fallback: all
    end
    imp_idx  = find(imp_mask);
    pf_sizes = arrayfun(@(r) size(r.F,1), results(imp_mask));
    [~, best] = max(pf_sizes);
    ri = imp_idx(best);
    fprintf('Auto-selected: results(%d)  Config=%s  [%s]  PF=%d pts\n', ...
        ri, results(ri).config_id, results(ri).config_desc, size(results(ri).F,1));
end

if ri < 1 || ri > n_res
    error('ri=%d out of range [1, %d]', ri, n_res);
end

r = results(ri);

if isempty(r.X) || size(r.X,1) == 0
    error('results(%d) has no Pareto X data. Re-run E1 after this fix.', ri);
end
if isempty(r.F) || size(r.F,1) == 0
    error('results(%d) has no Pareto F data.', ri);
end

n_pts = size(r.F, 1);

%% --- Auto-select Pareto point ---
if nargin < 3 || isempty(pi_idx)
    pi_idx = max(1, round(n_pts / 2));
end
pi_idx = min(max(1, pi_idx), n_pts);

%% --- Print Pareto front table ---
fprintf('\nPareto front for results(%d) — Config %s (%s):\n', ...
    ri, r.config_id, r.config_desc);
fprintf('  %5s  %8s  %8s\n', 'Idx', 'T(s)', 'E(kWh)');
for k = 1:n_pts
    sel = '';
    if k == pi_idx, sel = '  <-- selected'; end
    fprintf('  %5d  %8.2f  %8.4f%s\n', k, r.F(k,1), r.F(k,2), sel);
end
fprintf('\n');

X_sol = r.X(pi_idx, :);

%% --- Setup globals ---
route_file = r.route_file;
rs_file    = r.rs_file;
use_imp    = r.use_improved;

routePath = which(route_file);
rsPath    = which(rs_file);

if isempty(routePath)
    error(['Route file not found in path: %s\n' ...
        'Make sure you ran addpath(genpath(...)) for the project root.'], route_file);
end
if isempty(rsPath)
    error('Rollingstock file not found in path: %s', rs_file);
end

global use_improved driving_strategy  %#ok<NUSED>
use_improved     = use_imp;           %#ok<NASGU>
driving_strategy = "CC_CR";           %#ok<NASGU>
init_worker_globals(routePath, rsPath, struct(), use_imp);

%% --- Simulate ---
fprintf('Simulating Config %s, Pareto point %d/%d (T=%.1fs, E=%.4f kWh)...\n', ...
    r.config_id, pi_idx, n_pts, r.F(pi_idx,1), r.F(pi_idx,2));

if use_imp
    [running_inter, Total_E, s_out, v_out, vlim_out] = ...
        simulation_fun_CC_CR_improved(X_sol);
else
    [running_inter, Total_E, s_out, v_out, vlim_out] = ...
        simulation_fun_CC_CR_base(X_sol);
end

T_total = sum(running_inter);
fprintf('  Simulation result: T=%.2f s | E=%.4f kWh\n', T_total, Total_E);

%% --- Plot ---
fig_title = sprintf('Config %s — %s', r.config_id, r.config_desc);
figure('Name', fig_title, 'NumberTitle','off');
clf; hold on;

plot(s_out, v_out   * 3.6, 'b-',  'LineWidth', 2,   'DisplayName', 'Train speed');
plot(s_out, vlim_out* 3.6, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Speed limit');

xlabel('Distance (m)');
ylabel('Speed (km/h)');
title(sprintf('%s\nT = %.1f s  |  E = %.4f kWh  |  Pareto pt %d/%d', ...
    fig_title, T_total, Total_E, pi_idx, n_pts));
legend('Location','best');
grid on; box on;

end

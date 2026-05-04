%% run_E5_llm_advisor.m — Experiment E5: LLM Advisor Study
%
% Compares baseline configs against improved-strategy variants that use a
% file-based LLM advisor to generate route-aware smart-init priors.
%
% Variants:
%   A      — original CC_CR + NSGA-II vanilla
%   B      — original CC_CR + BRL-SDE-NSGA-II
%   C      — improved CC_CR + NSGA-II vanilla
%   D      — improved CC_CR + BRL-SDE-NSGA-II
%   C_LLM  — C + route-aware LLM advisor priors
%   D_LLM  — D + route-aware LLM advisor priors
%
% Workflow:
%   1. Export section-level route summary + prompt for an external LLM
%   2. Read JSON response file with advisor priors
%   3. Run the 6 configs with identical optimisation budget
%   4. Save comparable HV/IGD/feasibility/energy metrics

clc;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% ===== SETTINGS =====
if exist('ACTIVE_SEG','var') && ~isempty(ACTIVE_SEG)
    ROUTE_FILE = ACTIVE_SEG.file;
    T_TARGET   = ACTIVE_SEG.T_sched;
    N_RUNS     = ACTIVE_NRUNS;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', ACTIVE_SEG.name);
    SEG_LABEL  = ACTIVE_SEG.name;
else
    ROUTE_FILE = 'Guangzhou_Line7_IS02_1.120-3.028km.mat';
    T_TARGET   = 170;
    N_RUNS     = 10;
    OUT_DIR    = fullfile(fileparts(mfilename('fullpath')), '..', 'experiment_results', 'IS02');
    SEG_LABEL  = 'IS02';
end

RS_FILE    = 'rollingstock_Guangzhou_L7.m';
POP_SIZE   = 300;
ITERATIONS = 300;
N_WORKERS  = 4;
SEARCH_TLIM = T_TARGET * 1.50;
AUTO_GENERATE_RESPONSE = false;  % false = stop after prompt, then place response JSON and rerun

REQUEST_FILE  = fullfile(OUT_DIR, 'E5_llm_request.json');
PROMPT_FILE   = fullfile(OUT_DIR, 'E5_llm_prompt.txt');
RESPONSE_FILE = fullfile(OUT_DIR, 'E5_llm_response.json');

CONFIGS = struct( ...
    'id',             {'A', 'B', 'C', 'D', 'C_LLM', 'D_LLM'}, ...
    'use_improved',   {false, false, true, true, true, true}, ...
    'nsga2_variant',  {'original', 'rl_sde', 'original', 'rl_sde', 'original', 'rl_sde'}, ...
    'use_llm_advisor',{false, false, false, false, true, true}, ...
    'desc',           {'Original CC_CR + Original NSGA-II', ...
                       'Original CC_CR + BRL-SDE-NSGA-II', ...
                       'Improved CC_CR + Original NSGA-II', ...
                       'Improved CC_CR + BRL-SDE-NSGA-II', ...
                       'Improved CC_CR + Original NSGA-II + LLM advisor', ...
                       'Improved CC_CR + BRL-SDE-NSGA-II + LLM advisor'} );

if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

%% ===== GLOBALS + SETUP =====
global use_improved nsga2_variant pop_size iterations dimension var
global vel_profile time_obj_max parallel_use show_progress driving_strategy

show_progress    = false;
parallel_use     = true;
driving_strategy = "CC_CR";
time_obj_max     = SEARCH_TLIM;

cfg = struct('route_file', ROUTE_FILE, 'rollingstock_file', RS_FILE, ...
    'driving_strategy', "CC_CR", 'pop_size', POP_SIZE, 'iterations', ITERATIONS, ...
    'use_improved', false, 'nsga2_variant', 'original', 'parallel_use', true, 'sim', struct());
info = setup_project(cfg);

% Build route summary from improved sectioning before starting the pool.
init_worker_globals(info.routePath, info.rollingstockPath, cfg.sim, true);
advisor_request = build_llm_advisor_request(SEG_LABEL, T_TARGET, SEARCH_TLIM);
write_llm_advisor_files(advisor_request, REQUEST_FILE, PROMPT_FILE, RESPONSE_FILE);

if exist(RESPONSE_FILE, 'file') ~= 2
    if AUTO_GENERATE_RESPONSE
        auto_response = build_auto_llm_advisor_response(advisor_request);
        write_llm_advisor_response(auto_response, RESPONSE_FILE);
        fprintf('\n[E5] No advisor response found. Auto-generated route-aware response.\n');
        fprintf('  Request JSON : %s\n', REQUEST_FILE);
        fprintf('  Prompt file  : %s\n', PROMPT_FILE);
        fprintf('  Response JSON: %s\n\n', RESPONSE_FILE);
    else
        fprintf('\n[E5] LLM advisor request prepared.\n');
        fprintf('  Request JSON : %s\n', REQUEST_FILE);
        fprintf('  Prompt file  : %s\n', PROMPT_FILE);
        fprintf('  Response JSON: %s\n', RESPONSE_FILE);
        fprintf('Add the advisor response JSON, then rerun run_E5_llm_advisor.\n\n');
        return;
    end
end

advisor_profile = load_llm_advisor_profile(RESPONSE_FILE, var);

[~, E_BASELINE] = flat_out_baseline_rk4(info.routePath, info.rollingstockPath);
fprintf('E_baseline=%.4f kWh | T_target=%d s | Segment=%s\n', E_BASELINE, T_TARGET, SEG_LABEL);

setup_parallel_pool(cfg, info);
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', N_WORKERS);
    pool = gcp('nocreate');
end

%% ===== RUN ALL CONFIGS =====
n_cfg        = numel(CONFIGS);
pops_all     = cell(n_cfg, N_RUNS);
t_mat        = nan(n_cfg, N_RUNS);
dim_all      = zeros(n_cfg, 1);
union_PF     = cell(n_cfg, 1);

fprintf('\n===== E5 Run Plan =====\n');
for ci = 1:n_cfg
    fprintf('  %-6s | pop=%d | gen=%d | runs=%d | advisor=%d\n', ...
        CONFIGS(ci).id, POP_SIZE, ITERATIONS, N_RUNS, CONFIGS(ci).use_llm_advisor);
end
fprintf('=======================\n\n');

for ci = 1:n_cfg
    Cfg = CONFIGS(ci);
    fprintf('\n=== Config %s: %s ===\n', Cfg.id, Cfg.desc);

    use_improved  = Cfg.use_improved;
    nsga2_variant = Cfg.nsga2_variant;

    init_worker_globals(info.routePath, info.rollingstockPath, cfg.sim, Cfg.use_improved);
    wait(parfevalOnAll(@init_worker_globals, 0, info.routePath, info.rollingstockPath, cfg.sim, Cfg.use_improved));

    dim_all(ci) = dimension;
    pop_size    = POP_SIZE;
    iterations  = ITERATIONS;

    wait(parfevalOnAll(@push_worker_globals, 0, POP_SIZE, ITERATIONS, SEARCH_TLIM, Cfg.nsga2_variant));

    if Cfg.use_llm_advisor
        push_llm_advisor_globals(true, advisor_profile);
        wait(parfevalOnAll(@push_llm_advisor_globals, 0, true, advisor_profile));
    else
        push_llm_advisor_globals(false, struct());
        wait(parfevalOnAll(@push_llm_advisor_globals, 0, false, struct()));
    end

    vp_local = vel_profile;
    pops_ci  = cell(N_RUNS,1);
    times_ci = zeros(N_RUNS,1);

    parfor (run_id = 1:N_RUNS, N_WORKERS)
        rng(run_id, 'twister');
        t0_ = tic;
        [p,~,~] = nsga2_main(vp_local);
        times_ci(run_id) = toc(t0_);
        pops_ci{run_id}  = p;
    end

    F_union = zeros(0,2);
    for run_id = 1:N_RUNS
        pops_all{ci, run_id} = pops_ci{run_id};
        t_mat(ci, run_id)    = times_ci(run_id);

        pop_r = pops_ci{run_id};
        if isempty(pop_r) || size(pop_r,2) < dim_all(ci)+2
            continue;
        end
        F = double(pop_r(:, dim_all(ci)+1:dim_all(ci)+2));
        valid = F(:,1) < 1e5 & F(:,2) < 1e5;
        if any(valid)
            F_union = [F_union; nd_filter(F(valid,:))]; %#ok<AGROW>
        end
    end
    union_PF{ci} = nd_filter(F_union);
    fprintf('  Done | dim=%d | union PF=%d pts | mean runtime=%.1fs\n', ...
        dim_all(ci), size(union_PF{ci},1), mean(times_ci));
end

%% ===== GLOBAL REFERENCE FRONT =====
allFronts = zeros(0,2);
for ci = 1:n_cfg
    if ~isempty(union_PF{ci})
        allFronts = [allFronts; union_PF{ci}]; %#ok<AGROW>
    end
end
F_ref_global = nd_filter(allFronts);
fprintf('\nGlobal reference front: %d pts\n', size(F_ref_global,1));

%% ===== METRICS + RESULT TABLE =====
hv_mat  = nan(n_cfg, N_RUNS);
igd_mat = nan(n_cfg, N_RUNS);
fr_mat  = nan(n_cfg, N_RUNS);
sv_mat  = nan(n_cfg, N_RUNS);
sp_mat  = nan(n_cfg, N_RUNS);
E_mat   = nan(n_cfg, N_RUNS);
rows    = {};

tmpl_r = struct('solver','','seed',0,'pop_size',0,'dim',0,'iterations',0, ...
    'runtime',NaN,'AllX',zeros(0,0),'AllF',zeros(0,2),'X',zeros(0,0),'F',zeros(0,2), ...
    'HV',NaN,'IGD',NaN,'Spread',NaN,'Nf1',0,'nsga2_variant','','config_id','', ...
    'config_desc','','use_improved',false,'use_llm_advisor',false, ...
    'route_file','','rs_file','');
results = repmat(tmpl_r, 0, 1);

for ci = 1:n_cfg
    d = dim_all(ci);
    for run_id = 1:N_RUNS
        pop_r = pops_all{ci, run_id};
        if isempty(pop_r) || size(pop_r,2) < d+2
            rows{end+1} = {CONFIGS(ci).id, CONFIGS(ci).use_llm_advisor, run_id, NaN, NaN, NaN, NaN, NaN, t_mat(ci,run_id)}; %#ok<AGROW>
            continue;
        end

        F_all = double(pop_r(:, d+1:d+2));
        [hv_mat(ci,run_id), igd_mat(ci,run_id), sp_mat(ci,run_id), fr_mat(ci,run_id)] = ...
            compute_metrics(F_all, F_ref_global, T_TARGET);

        valid = F_all(:,1) < 1e5 & F_all(:,2) < 1e5;
        F_nd  = nd_filter(F_all(valid,:));

        E_best = NaN;
        if ~isempty(F_nd)
            feas = F_nd(F_nd(:,1) <= T_TARGET, :);
            if ~isempty(feas)
                E_best = min(feas(:,2));
            end
        end
        E_mat(ci, run_id)  = E_best;
        sv_mat(ci, run_id) = energy_saving_percent(E_best, E_BASELINE);

        rows{end+1} = {CONFIGS(ci).id, CONFIGS(ci).use_llm_advisor, run_id, ...
            hv_mat(ci,run_id), igd_mat(ci,run_id), fr_mat(ci,run_id), E_best, sv_mat(ci,run_id), t_mat(ci,run_id)}; %#ok<AGROW>

        r = tmpl_r;
        r.solver          = 'nsga2';
        r.seed            = run_id;
        r.pop_size        = POP_SIZE;
        r.dim             = d;
        r.iterations      = ITERATIONS;
        r.runtime         = t_mat(ci, run_id);
        r.HV              = hv_mat(ci, run_id);
        r.IGD             = igd_mat(ci, run_id);
        r.Spread          = sp_mat(ci, run_id);
        r.nsga2_variant   = CONFIGS(ci).nsga2_variant;
        r.config_id       = CONFIGS(ci).id;
        r.config_desc     = CONFIGS(ci).desc;
        r.use_improved    = CONFIGS(ci).use_improved;
        r.use_llm_advisor = CONFIGS(ci).use_llm_advisor;
        r.route_file      = ROUTE_FILE;
        r.rs_file         = RS_FILE;

        r.AllX = double(pop_r(:, 1:d));
        r.AllF = F_all;
        if any(valid)
            X_valid = r.AllX(valid, :);
            [r.F, nd_idx] = nd_filter_idx(F_all(valid,:));
            r.X = X_valid(nd_idx, :);
            r.Nf1 = size(r.F, 1);
        end
        results(end+1,1) = r; %#ok<AGROW>
    end
end

%% ===== SAVE OUTPUTS =====
T_csv = cell2table(vertcat(rows{:}), ...
    'VariableNames', {'config','use_llm_advisor','run_id','hv','igd','feasible','E_best_kWh','energy_saving_pct','runtime_s'});
writetable(T_csv, fullfile(OUT_DIR, 'E5_results.csv'));
meta = struct( ...
    'segment_name', SEG_LABEL, ...
    'T_target_s', T_TARGET, ...
    'E_baseline_kWh', E_BASELINE, ...
    'pop_size', POP_SIZE, ...
    'iterations', ITERATIONS, ...
    'n_runs', N_RUNS, ...
    'n_workers', N_WORKERS, ...
    'search_tlim_s', SEARCH_TLIM, ...
    'route_file', ROUTE_FILE, ...
    'rs_file', RS_FILE, ...
    'response_file', RESPONSE_FILE);
metrics_table = T_csv;
save(fullfile(OUT_DIR, 'E5_all_configs.mat'), 'results', 'advisor_profile', 'advisor_request', 'meta', 'metrics_table');

%% ===== SUMMARY =====
fid = fopen(fullfile(OUT_DIR, 'E5_summary.txt'), 'w');
fprintf(fid, 'Experiment E5 — LLM Advisor Study\n');
fprintf(fid, 'Segment: %s | T_target=%d s | E_baseline=%.4f kWh\n', SEG_LABEL, T_TARGET, E_BASELINE);
fprintf(fid, 'Pop=%d | Gen=%d | Runs=%d | Workers=%d\n', POP_SIZE, ITERATIONS, N_RUNS, N_WORKERS);
fprintf(fid, 'Response file: %s\n\n', RESPONSE_FILE);
fprintf(fid, '%-12s  %8s  %8s  %8s  %10s  %10s\n', ...
    'Config', 'HV_mean', 'HV_std', 'IGD_mean', 'Feas_mean', 'EnSav_pct');

for ci = 1:n_cfg
    line = sprintf('%-12s  %8.4f  %8.4f  %8.4f  %10.3f  %10.2f', ...
        CONFIGS(ci).id, mean(hv_mat(ci,:),'omitnan'), std(hv_mat(ci,:),'omitnan'), ...
        mean(igd_mat(ci,:),'omitnan'), mean(fr_mat(ci,:),'omitnan'), mean(sv_mat(ci,:),'omitnan'));
    disp(line);
    fprintf(fid, '%s\n', line);
end

fprintf(fid, '\nWilcoxon signed-rank test (HV, paired seeds) — D_LLM vs others:\n');
ref_idx = find(strcmp({CONFIGS.id}, 'D_LLM'), 1, 'first');
for ci = 1:n_cfg
    if ci == ref_idx
        continue;
    end
    x = hv_mat(ref_idx, :)';
    y = hv_mat(ci, :)';
    ok = ~isnan(x) & ~isnan(y);
    if sum(ok) >= 5 && exist('signrank','file') == 2
        p_val = signrank(x(ok), y(ok));
    else
        p_val = wilcoxon_srtest(x(ok), y(ok));
    end
    fprintf(fid, '  D_LLM vs %-6s : p = %.4f%s\n', CONFIGS(ci).id, p_val, ternary(p_val < 0.05, '  *', ''));
end
fclose(fid);

figure(5101); clf;
hv_means = mean(hv_mat, 2, 'omitnan');
hv_stds  = std(hv_mat, 0, 2, 'omitnan');
b = bar(hv_means, 0.65, 'FaceColor', 'flat');
hold on;
errorbar(1:n_cfg, hv_means, hv_stds, 'k.', 'LineWidth', 1.5, 'CapSize', 8);
b.CData = [0.55 0.55 0.55; 0.20 0.45 0.85; 0.95 0.50 0.05; 0.08 0.78 0.18; 0.95 0.70 0.25; 0.15 0.55 0.15];
set(gca, 'XTick', 1:n_cfg, 'XTickLabel', {CONFIGS.id}, 'FontSize', 10);
xlabel('Configuration');
ylabel('Hypervolume (HV)');
title(sprintf('E5: LLM Advisor Comparison — %s', SEG_LABEL));
grid on; box on;
saveas(gcf, fullfile(OUT_DIR, 'E5_hv_bar.png'));

fprintf('\nE5 complete → %s\n', OUT_DIR);

%% ===== LOCAL HELPERS =====
function request = build_llm_advisor_request(seg_label, t_target, search_tlim)
    global vel_profile gradient var station_info dimension

    nSec = size(vel_profile, 1) - 1;
    sections = repmat(struct( ...
        'section_id', 0, 'start_km', 0, 'end_km', 0, 'length_m', 0, ...
        'speed_limit_kmh', 0, 'policy_type', '', 'expected_dv', 0, ...
        'mean_gradient_permille', 0, 'min_gradient_permille', 0, ...
        'max_gradient_permille', 0, 'distance_to_next_station_m', NaN), nSec, 1);

    g_bp    = gradient(:,1);
    g_val   = gradient(:,2);
    g_start = g_bp(1:end-1);
    g_end   = g_bp(2:end);
    g_seg   = g_val(1:end-1);

    for i = 1:nSec
        a = vel_profile(i,1);
        b = vel_profile(i+1,1);
        overlap = max(0, min(g_end, b) - max(g_start, a));
        idx = overlap > 0;
        if any(idx)
            w = overlap(idx);
            g_here = g_seg(idx);
            g_mean = sum(g_here .* w) / max(sum(w), eps);
            g_min  = min(g_here);
            g_max  = max(g_here);
        else
            g_mean = interp1(gradient(:,1), gradient(:,2), (a+b)/2, 'linear', 'extrap');
            g_min  = g_mean;
            g_max  = g_mean;
        end

        next_station = station_info(find(station_info(:,1) > b, 1, 'first'), 1); %#ok<FNDSB>
        if isempty(next_station)
            dist_station = NaN;
        else
            dist_station = max(0, (next_station - b) * 1000);
        end

        sections(i).section_id = i;
        sections(i).start_km = a;
        sections(i).end_km = b;
        sections(i).length_m = (b - a) * 1000;
        sections(i).speed_limit_kmh = vel_profile(i,2);
        if var(i) == 3
            sections(i).policy_type = 'downhill';
            sections(i).expected_dv = 2;
        else
            sections(i).policy_type = 'normal';
            sections(i).expected_dv = min(var(i), 2);
        end
        sections(i).mean_gradient_permille = g_mean;
        sections(i).min_gradient_permille = g_min;
        sections(i).max_gradient_permille = g_max;
        sections(i).distance_to_next_station_m = dist_station;
    end

    request = struct();
    request.segment_name = seg_label;
    request.t_target_s = t_target;
    request.search_tlim_s = search_tlim;
    request.dimension = dimension;
    request.n_sections = nSec;
    request.strategy_notes = ['Improved CC_CR uses section-level priors. ' ...
        'For policy_type="normal", output cruise_ratio and coast_ratio. ' ...
        'For policy_type="downhill", output high_ratio and low1_ratio. ' ...
        'Do not output MID; it is already computed adaptively inside the simulator.'];
    request.sections = sections;
end

function write_llm_advisor_files(request, request_file, prompt_file, response_file)
    req_text = jsonencode(request);
    fid = fopen(request_file, 'w');
    fprintf(fid, '%s', req_text);
    fclose(fid);

    prompt = sprintf([ ...
        'You are designing route-aware priors for an improved CC_CR train-driving strategy.\n\n' ...
        'Task:\n' ...
        'Return JSON only, with exactly one top-level key: "sections".\n' ...
        'Each entry must match a section_id from the route summary.\n' ...
        'Rules:\n' ...
        '- For policy_type="normal": return cruise_ratio and coast_ratio.\n' ...
        '- For policy_type="downhill": return high_ratio and low1_ratio.\n' ...
        '- All ratios must be in (0,1].\n' ...
        '- coast_ratio must not imply coasting above cruising.\n' ...
        '- low1_ratio must not imply LOW1 above HIGH.\n' ...
        '- MID is already gradient-adaptive in the simulator, so do not emit MID.\n' ...
        '- Optimize for meeting T_target while reducing energy, especially on downhill sections.\n\n' ...
        'Expected JSON schema:\n' ...
        '{"sections":[{"section_id":1,"policy_type":"normal","cruise_ratio":0.92,"coast_ratio":0.74},' ...
        '{"section_id":2,"policy_type":"downhill","high_ratio":0.85,"low1_ratio":0.67}]}\n\n' ...
        'Route summary JSON:\n%s\n\n' ...
        'Write your response JSON to:\n%s\n'], req_text, response_file);

    fid = fopen(prompt_file, 'w');
    fprintf(fid, '%s', prompt);
    fclose(fid);
end

function response = build_auto_llm_advisor_response(request)
    sections = request.sections;
    nSec = numel(sections);
    response = struct();
    base_entry = struct( ...
        'section_id', 0, ...
        'policy_type', '', ...
        'cruise_ratio', NaN, ...
        'coast_ratio', NaN, ...
        'high_ratio', NaN, ...
        'low1_ratio', NaN);
    response.sections = repmat(base_entry, nSec, 1);

    for i = 1:nSec
        sec = sections(i);
        dist_to_station = sec.distance_to_next_station_m;
        if isempty(dist_to_station) || ~isfinite(dist_to_station)
            dist_to_station = sec.length_m;
        end

        long_factor = min(max(dist_to_station / 1500, 0), 1);
        uphill_factor = min(max(sec.mean_gradient_permille, 0) / 15, 1);
        downhill_factor = min(max(-sec.mean_gradient_permille, 0) / 25, 1);

        out = struct();
        out = base_entry;
        out.section_id = sec.section_id;
        out.policy_type = sec.policy_type;

        if strcmpi(sec.policy_type, 'downhill')
            high_ratio = 0.82 - 0.08 * downhill_factor + 0.03 * long_factor;
            low1_ratio = 0.62 - 0.10 * downhill_factor + 0.04 * long_factor;
            high_ratio = clamp_ratio(high_ratio);
            low1_ratio = min(clamp_ratio(low1_ratio), high_ratio);
            out.high_ratio = high_ratio;
            out.low1_ratio = low1_ratio;
        else
            cruise_ratio = 0.84 + 0.07 * long_factor + 0.04 * uphill_factor - 0.03 * downhill_factor;
            coast_ratio  = 0.68 + 0.06 * long_factor - 0.03 * uphill_factor + 0.02 * downhill_factor;
            cruise_ratio = clamp_ratio(cruise_ratio);
            coast_ratio  = min(clamp_ratio(coast_ratio), cruise_ratio);
            out.cruise_ratio = cruise_ratio;
            out.coast_ratio  = coast_ratio;
        end

        response.sections(i) = out;
    end
end

function write_llm_advisor_response(response, file_path)
    fid = fopen(file_path, 'w');
    fprintf(fid, '%s', jsonencode(response));
    fclose(fid);
end

function profile = load_llm_advisor_profile(file_path, var)
    raw = fileread(file_path);
    data = jsondecode(raw);
    if ~isstruct(data) || ~isfield(data, 'sections')
        error('LLM advisor response must contain top-level key: sections');
    end

    sections = data.sections;
    if ~isstruct(sections)
        error('LLM advisor sections must decode to a struct array.');
    end
    sections = sections(:);
    if numel(sections) ~= numel(var)
        error('Advisor section count mismatch: expected %d, got %d', numel(var), numel(sections));
    end

    for i = 1:numel(var)
        if ~isfield(sections(i), 'section_id') || double(sections(i).section_id) ~= i
            error('Advisor section_id mismatch at entry %d', i);
        end

        switch var(i)
            case 3
                assert(isfield(sections(i), 'high_ratio') && isfield(sections(i), 'low1_ratio'), ...
                    'Downhill section %d must define high_ratio and low1_ratio', i);
                sections(i).high_ratio = clamp_ratio(sections(i).high_ratio);
                sections(i).low1_ratio = min(clamp_ratio(sections(i).low1_ratio), sections(i).high_ratio);
                sections(i).policy_type = 'downhill';

            otherwise
                assert(isfield(sections(i), 'cruise_ratio') && isfield(sections(i), 'coast_ratio'), ...
                    'Normal section %d must define cruise_ratio and coast_ratio', i);
                sections(i).cruise_ratio = clamp_ratio(sections(i).cruise_ratio);
                sections(i).coast_ratio  = min(clamp_ratio(sections(i).coast_ratio), sections(i).cruise_ratio);
                sections(i).policy_type = 'normal';
        end
    end

    profile = struct('sections', sections, 'source_file', file_path);
end

function val = clamp_ratio(val)
    val = min(max(double(val), 0.05), 0.99);
end

function F_nd = nd_filter(F)
    if isempty(F)
        F_nd = zeros(0,2);
        return;
    end
    keep = true(size(F,1),1);
    for i = 1:size(F,1)
        if ~keep(i)
            continue;
        end
        for j = 1:size(F,1)
            if i == j || ~keep(j)
                continue;
            end
            if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                keep(i) = false;
                break;
            end
        end
    end
    F_nd = F(keep,:);
    [~, ord] = sort(F_nd(:,1), 'ascend');
    F_nd = F_nd(ord,:);
end

function [F_nd, idx_nd] = nd_filter_idx(F)
    if isempty(F)
        F_nd = zeros(0,2);
        idx_nd = zeros(0,1);
        return;
    end
    keep = true(size(F,1),1);
    for i = 1:size(F,1)
        if ~keep(i)
            continue;
        end
        for j = 1:size(F,1)
            if i == j || ~keep(j)
                continue;
            end
            if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                keep(i) = false;
                break;
            end
        end
    end
    idx_nd = find(keep);
    F_nd = F(idx_nd,:);
    [~, ord] = sort(F_nd(:,1), 'ascend');
    F_nd = F_nd(ord,:);
    idx_nd = idx_nd(ord);
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
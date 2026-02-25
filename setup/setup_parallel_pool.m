function setup_parallel_pool(cfg, info)
%SETUP_PARALLEL_POOL Start/resize pool, attach files, init globals on workers.
%
% Required:
%   init_worker_globals(routePath, rollingstockPath, sim_in, use_improved_in)
%
% Notes:
% - We attach ONLY critical files to workers for stability.
% - When benchmark switches strategy, run_benchmark_compare() will call
%   init_worker_globals again (broadcast).

    project_root = info.project_root;

    % ---------- Decide worker count ----------
    c = parcluster('local');
    nWorkers = c.NumWorkers;
    if isfield(cfg,'leave_one_core') && cfg.leave_one_core
        nWorkers = max(1, nWorkers - 1);
    end

    % ---------- Start/resize pool ----------
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers ~= nWorkers
        delete(pool);
        pool = [];
    end
    if isempty(pool)
        parpool('local', nWorkers);
    end
    pool = gcp;
    fprintf('Pool aktif: %d workers | Profile: %s\n', pool.NumWorkers, pool.Cluster.Profile);

    % ---------- Attach essential files ----------
    filesToAttach = {};

    % route + rollingstock files
    if isfield(info,'routePath') && ~isempty(info.routePath), filesToAttach{end+1} = info.routePath; end
    if isfield(info,'rollingstockPath') && ~isempty(info.rollingstockPath), filesToAttach{end+1} = info.rollingstockPath; end

    % simulators (CC_CR only here; add CC/CR if you use them)
    simFiles = {'simulation_fun_CC_CR.m','simulation_fun_CC_CR_base.m','simulation_fun_CC_CR_improved.m'};
    for k = 1:numel(simFiles)
        p = which(simFiles{k});
        if ~isempty(p)
            filesToAttach{end+1} = p; %#ok<AGROW>
        else
            error('Simulator file not found on path: %s', simFiles{k});
        end
    end

    % decision var dispatcher + variants + worker init
    mustFiles = {'calculate_pop.m', 'decision_var_NO.m','decision_var_NO_base.m','decision_var_NO_improved.m', 'init_worker_globals.m'};
    for k = 1:numel(mustFiles)
        p = which(mustFiles{k});
        if isempty(p)
            error('Required file not found on path: %s', mustFiles{k});
        end
        filesToAttach{end+1} = p; %#ok<AGROW>
    end

    % solver files (optional)
    maybeSolvers = {'nsga2_main.m','mopso_main.m','spea2_main.m','moead_main.m'};
    for k = 1:numel(maybeSolvers)
        p = which(maybeSolvers{k});
        if ~isempty(p)
            filesToAttach{end+1} = p; %#ok<AGROW>
        end
    end

    filesToAttach = unique(filesToAttach, 'stable');

    addAttachedFiles(pool, filesToAttach);

    % ---------- Ensure project path on workers ----------
    fPath = parfevalOnAll(@() addpath(genpath(project_root)), 0);
    wait(fPath);

    % ---------- Initialize globals on workers ----------
    if ~isfield(cfg,'sim'), cfg.sim = struct(); end
    if ~isfield(cfg,'use_improved'), cfg.use_improved = false; end

    fInit = parfevalOnAll(@init_worker_globals, 0, info.routePath, info.rollingstockPath, cfg.sim, cfg.use_improved);
    try
        fetchOutputs(fInit);
    catch ME
        fprintf('\n======================================================\n');
        fprintf('!!! FATAL ERROR SAAT INISIALISASI WORKER !!!\n');
        fprintf('======================================================\n');
        disp(ME.getReport());
        error('Worker gagal inisialisasi. Cek pesan error di atas bos!');
    end

    % ---------- Print PIDs (debug) ----------
    try
        fPID = parfevalOnAll(@() feature('getpid'), 1);
        pids = fetchOutputs(fPID);
        disp(pids);
    catch
        % ignore
    end
end

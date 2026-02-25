function setup_globals_parallel(routeMatFile, rollingstockPath, use_improved_in)
%SETUP_GLOBALS_PARALLEL Legacy helper (optional).
% Prefer using init_worker_globals + setup_parallel_pool.
%
% This keeps behavior consistent with base vs improved selection.

    global vel_profile station_info gradient terminal
    global Mass lambda inertial_mass Davis gravity
    global max_speed
    global var dimension
    global use_improved

    global vel_profile_raw gradient_raw

    if nargin < 3 || isempty(use_improved_in)
        use_improved = false;
    else
        use_improved = logical(use_improved_in);
    end

    if isempty(gravity), gravity = 9.81; end

    % clear caches
    vel_profile_raw = [];
    gradient_raw    = [];

    % load route
    S = load(routeMatFile);
    vel_profile  = S.vel_profile;
    station_info = S.station_info;
    gradient     = S.gradient;
    terminal     = S.terminal;

    % rollingstock
    run(rollingstockPath);
    if isempty(max_speed)
        error('max_speed masih kosong. Cek rollingstock file.');
    end
    max_speed = max_speed(1);

    % decision vars
    var = zeros(1, length(vel_profile)-1);
    decision_var_NO;
    dimension = sum(var);
end

function init_worker_globals(routePath, rollingstockPath, sim_in, use_improved_in)
%INIT_WORKER_GLOBALS Initialize all globals on each worker.
% This keeps DIMENSION consistent across workers when switching base vs improved.

    global vel_profile station_info gradient terminal
    global Mass lambda inertial_mass Davis gravity
    global V1_traction V2_traction V3_traction
    global V1_brake V2_brake V3_brake
    global Max_tractive_power Max_brake_power
    global max_speed co_fric max_traction max_brake
    global max_accel_trac max_accel_brake

    global var dimension driving_strategy
    global use_improved

    % Raw caches used by improved decision_var dispatcher
    global vel_profile_raw gradient_raw

    % ---------- mode flag (passed from main) ----------
    if nargin < 4 || isempty(use_improved_in)
        use_improved = false;
    else
        use_improved = logical(use_improved_in);
    end

    % ---------- optional simopts (reserved) ----------
    %#ok<NASGU>
    if nargin < 3
        sim_in = struct();
    end

    % ---------- driving strategy ----------
    if isempty(driving_strategy)
        driving_strategy = "CC_CR";
    else
        driving_strategy = string(driving_strategy);
    end

    if isempty(gravity), gravity = 9.81; end

    % IMPORTANT: clear raw caches so each init is deterministic
    vel_profile_raw = [];
    gradient_raw    = [];

    % ---------- load route ----------
    S = load(routePath);
    if isfield(S,'vel_profile'),  vel_profile  = S.vel_profile;  end
    if isfield(S,'station_info'), station_info = S.station_info; end
    if isfield(S,'gradient'),     gradient     = S.gradient;     end
    if isfield(S,'terminal'),     terminal     = S.terminal;     end
    if isempty(vel_profile)
        error('[WORKER INIT] vel_profile missing after loading route: %s', routePath);
    end

    % ---------- load rollingstock ----------
    run(rollingstockPath);

    if exist('rotary','var') && ~isempty(rotary)
        lambda = rotary;
    end
    if isempty(lambda), lambda = 0.07; end
    lambda = lambda(1);

    if ~isempty(max_speed), max_speed = max_speed(1); end

    % ---------- sync var/dimension via dispatcher ----------
    var = zeros(1, length(vel_profile)-1);
    decision_var_NO;
    if logical(use_improved)
        % Gradient-adaptive MID: improved CC_CR downhill sections (var=3) use only 2 DVs.
        % MID is computed from gradient inside the simulator, not a decision variable.
        dimension = sum(min(var, 2));
    else
        dimension = sum(var);
    end

    fprintf('[WORKER INIT] use_improved=%d | strategy=%s | dim=%d\n', ...
        use_improved, upper(string(driving_strategy)), dimension);
end

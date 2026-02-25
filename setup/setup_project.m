function info = setup_project(cfg)
%SETUP_PROJECT Load route + rollingstock and compute decision dimension.
%
% This version supports TWO strategies (base vs improved) consistently:
%   - simulator selector: global use_improved
%   - decision variables: decision_var_NO.m dispatcher
%
% Required files in project root:
%   - simulation_fun_CC_CR.m          (WRAPPER)
%   - simulation_fun_CC_CR_base.m
%   - simulation_fun_CC_CR_improved.m
%   - decision_var_NO.m               (DISPATCHER)
%   - decision_var_NO_base.m
%   - decision_var_NO_improved.m

% ================= GLOBALS =================
global var pop_size iterations dimension
global vel_profile station_info gradient terminal
global Mass lambda inertial_mass Davis
global V1_traction V2_traction V3_traction
global V1_brake V2_brake V3_brake
global Max_tractive_power Max_brake_power
global max_speed co_fric gravity max_traction max_brake
global max_accel_trac max_accel_brake
global parallel_use driving_strategy

global use_improved

% ================= BASIC CFG =================
parallel_use     = logical(cfg.parallel_use);
driving_strategy = string(cfg.driving_strategy);

pop_size   = cfg.pop_size;
iterations = cfg.iterations;

if isfield(cfg,'use_improved')
    use_improved = logical(cfg.use_improved);
else
    use_improved = false;
end

% ================= LOCATE FILES =================
routePath = which(cfg.route_file);
if isempty(routePath)
    if exist(cfg.route_file,'file')==2
        routePath = fullfile(pwd, cfg.route_file);
    else
        error('Route file not found: %s', cfg.route_file);
    end
end

rsPath = which(cfg.rollingstock_file);
if isempty(rsPath)
    if exist(cfg.rollingstock_file,'file')==2
        rsPath = fullfile(pwd, cfg.rollingstock_file);
    else
        error('Rollingstock file not found: %s', cfg.rollingstock_file);
    end
end

% ================= LOAD ROUTE =================
S = load(routePath);
if isfield(S,'vel_profile'),  vel_profile  = S.vel_profile;  end
if isfield(S,'station_info'), station_info = S.station_info; end
if isfield(S,'gradient'),     gradient     = S.gradient;     end
if isfield(S,'terminal'),     terminal     = S.terminal;     end

if isempty(vel_profile)
    error('vel_profile not loaded from route file.');
end

% ================= LOAD ROLLINGSTOCK =================
run(rsPath);

% rotary -> lambda
if exist('rotary','var') && ~isempty(rotary)
    lambda = rotary;
end
if isempty(lambda)
    lambda = 0.07;
end
lambda = lambda(1);   % make scalar

% gravity
if isempty(gravity)
    gravity = 9.81;
end

% Ensure max_speed scalar if provided as vector
if ~isempty(max_speed)
    max_speed = max_speed(1);
end

% Derive traction/brake piecewise params if missing (keep your old logic)
if isempty(V1_traction) && exist('Power','var') && ~isempty(max_traction)
    V1_traction = Power/1000/max_traction;
end
if isempty(V2_traction) && ~isempty(max_speed)
    V2_traction = max_speed/3.6; % m/s
end
if isempty(V1_brake), V1_brake = V1_traction; end
if isempty(V2_brake), V2_brake = V2_traction; end

if isempty(Max_tractive_power) && exist('Power','var')
    Max_tractive_power = Power;
end
if isempty(Max_brake_power) && exist('Power','var')
    Max_brake_power = Power;
end

% ================= COMPUTE var/dimension =================
var = zeros(1, length(vel_profile)-1);

% decision_var_NO is a DISPATCHER (base vs improved) controlled by global use_improved
% It may refine vel_profile (improved mode).
decision_var_NO;

dimension = sum(var);

fprintf('Setup done | strategy=%s | use_improved=%d | dimension=%d | pop=%d | iter=%d\n', ...
    driving_strategy, use_improved, dimension, pop_size, iterations);

% ================= RETURN INFO =================
info = struct();
info.project_root     = pwd;
info.routePath        = routePath;
info.rollingstockPath = rsPath;
info.dimension        = dimension;
info.pop_size         = pop_size;
info.iterations       = iterations;
end

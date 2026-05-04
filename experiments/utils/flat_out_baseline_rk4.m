function [T_min, E_baseline] = flat_out_baseline_rk4(segment_file, rollingstock_file)
%FLAT_OUT_BASELINE_RK4  Compute minimum running time and baseline energy via RK4.
%
% Simulates the train from start to stop using:
%   - Phase 1: Full traction (max acceleration / max power)
%   - Phase 2: Full braking (max deceleration)
% Automatically switches from traction to braking at the optimal point
% so that the train decelerates to rest exactly at the terminal station.
%
% Physics integration: Runge-Kutta 4, dt = 0.1 s
%
% Inputs:
%   segment_file     : path or name of .mat file containing
%                      vel_profile, station_info, gradient, terminal
%   rollingstock_file: path or name of .m script defining rolling stock
%                      parameters (Mass, lambda, Davis, traction/brake data)
%
% Outputs:
%   T_min     : minimum running time (seconds)
%   E_baseline: traction energy consumption (kWh), accounting for regen braking

    % ===== LOAD SEGMENT DATA =====
    seg_path = resolve_path(segment_file);
    S = load(seg_path);
    vel_profile_s  = S.vel_profile;   % N×2 [dist_km, speed_limit_km/h]
    station_info_s = S.station_info;  % M×2 [station_dist_km, dwell_s]
    gradient_raw   = S.gradient;      % K×2 [dist_km, grade_permille]

    % Terminal station distance (km → m)
    S_end = max(station_info_s(:,1)) * 1000;   % m

    % ===== LOAD ROLLING STOCK =====
    rs_path = resolve_path(rollingstock_file);
    tmp = run_rs_script(rs_path);

    Mass_t            = tmp.Mass;
    lambda_r          = tmp.lambda;
    Davis_rs          = tmp.Davis;
    Max_tractive_power= tmp.Max_tractive_power / 1000;  % W → kW
    Max_brake_power   = tmp.Max_brake_power   / 1000;
    max_accel_trac_rs = tmp.max_accel_trac;
    max_accel_brake_rs= tmp.max_accel_brake;
    inertial_mass_t   = Mass_t * (1 + lambda_r);
    REGEN_EFF         = 0.6;

    % ===== BUILD INTERPOLANTS =====
    vl_s = vel_profile_s(:,1) * 1000;   % m
    vl_v = vel_profile_s(:,2) / 3.6;    % m/s
    vlim_fn = @(s) interp1(vl_s, vl_v, s, 'previous', vl_v(end));

    gr_s = gradient_raw(:,1) * 1000;
    gr_g = gradient_raw(:,2);
    if numel(gr_s) > 1
        grade_fn = @(s) interp1(gr_s, gr_g, s, 'previous', 0);
    else
        grade_fn = @(~) 0;
    end

    % ===== PRECOMPUTE BRAKING DISTANCE LUT =====
    v_max_ms  = max(vl_v);
    dv_lut    = 0.5;
    v_lut     = 0 : dv_lut : v_max_ms + 1;
    d_brk_lut = zeros(size(v_lut));
    DT_BRK    = 0.05;

    for vi = length(v_lut):-1:2
        v0 = v_lut(vi);
        v  = v0; d = 0;
        while v > 0.05
            a_brk = compute_brk_accel(v, max_accel_brake_rs, Max_brake_power, ...
                                       inertial_mass_t, Davis_rs, 0);
            v_new = max(0, v + a_brk * DT_BRK);
            ds    = (v + v_new) / 2 * DT_BRK;
            d = d + ds;
            v = v_new;
        end
        d_brk_lut(vi) = d;
    end
    brk_dist_fn = @(v) interp1(v_lut, d_brk_lut, min(v, v_max_ms+1), 'linear', 0);

    % ===== RK4 FORWARD SIMULATION =====
    DT   = 0.1;
    s    = 0.0;
    v    = 0.0;
    t    = 0.0;
    E_trac = 0.0;
    E_brk  = 0.0;
    phase  = 1;
    MAX_T  = 600;

    while s < S_end - 0.5 && t < MAX_T
        if phase == 1
            if (S_end - s) <= brk_dist_fn(v) * 1.02
                phase = 2;
            end
        end

        vlim_cur = vlim_fn(s);

        if phase == 1
            f_a = @(ss,vv) compute_trac_accel(vv, max_accel_trac_rs, ...
                Max_tractive_power, inertial_mass_t, Davis_rs, grade_fn(ss), vlim_fn(ss));
        else
            f_a = @(ss,vv) compute_brk_accel(vv, max_accel_brake_rs, ...
                Max_brake_power, inertial_mass_t, Davis_rs, grade_fn(ss));
        end

        [s_new, v_new] = rk4_step(s, v, DT, f_a);
        v_new = max(0, min(v_new, vlim_cur * 1.001));

        v_avg = (v + v_new) / 2;
        if v_avg > 0.1
            if phase == 1
                F_kn  = min(max_accel_trac_rs * inertial_mass_t, Max_tractive_power / v_avg);
                E_trac = E_trac + F_kn * v_avg * DT / 3600;
            else
                F_kn  = min(max_accel_brake_rs * inertial_mass_t, Max_brake_power / v_avg);
                E_brk  = E_brk  + F_kn * v_avg * DT / 3600;
            end
        end

        s = s_new; v = v_new; t = t + DT;

        if s >= S_end || (phase==2 && v < 0.1 && s > S_end*0.85)
            break;
        end
    end

    T_min      = t;
    E_baseline = max(0, E_trac - REGEN_EFF * E_brk);

    fprintf('[flat_out_rk4] %s | T_min=%.2f s | E_baseline=%.4f kWh\n', ...
        segment_name_from_path(segment_file), T_min, E_baseline);
end

% =========================================================================
%  PHYSICS HELPERS
% =========================================================================

function a = compute_trac_accel(v, a_trac_max, P_max_kw, m_t, Davis, grade_ppm, vlim)
    if v < 0.1
        a_trac = a_trac_max;
    else
        a_trac = min(a_trac_max, P_max_kw / (m_t * v));
    end
    a_resist = (Davis(1) + Davis(2)*v + Davis(3)*v^2) / m_t;
    a_grad   = 9.81 * grade_ppm / 1000;
    a = a_trac - a_resist + a_grad;
    if v >= vlim - 0.1
        a = min(0, a);
    end
end

function a = compute_brk_accel(v, a_brk_max, P_brk_kw, m_t, Davis, grade_ppm)
    if v < 0.1
        a_brk = -a_brk_max;
    else
        a_brk = -min(a_brk_max, P_brk_kw / (m_t * v));
    end
    a_resist = (Davis(1) + Davis(2)*v + Davis(3)*v^2) / m_t;
    a_grad   = 9.81 * grade_ppm / 1000;
    a = a_brk - a_resist + a_grad;
end

function [s_new, v_new] = rk4_step(s, v, dt, f_a)
    k1v = f_a(s, v);
    k2v = f_a(s + dt/2*v,              v + dt/2*k1v);
    k3v = f_a(s + dt/2*(v+dt/2*k1v),  v + dt/2*k2v);
    k4v = f_a(s + dt*(v+dt*k2v),       v + dt*k3v);
    v_new = max(0, v + dt/6*(k1v + 2*k2v + 2*k3v + k4v));
    s_new = s + dt/6*(v + 2*(v+dt/2*k1v) + 2*(v+dt/2*k2v) + (v+dt*k3v));
end

% =========================================================================
%  FILE / PATH HELPERS
% =========================================================================

function p = resolve_path(f)
    if exist(f,'file') == 2
        p = which(f);
        if isempty(p), p = f; end
        return;
    end
    p = which(f);
    if isempty(p)
        error('File not found: %s', f);
    end
end

function name = segment_name_from_path(f)
    [~, fname, ~] = fileparts(f);
    tok = regexp(fname, 'IS\d+', 'match');
    if ~isempty(tok)
        name = tok{1};
    else
        name = fname(1:min(20,numel(fname)));
    end
    name = strrep(name, '-', '_');
end

function vars = run_rs_script(rs_path)
    Mass = []; lambda = []; Davis = []; %#ok<NASGU>
    Max_tractive_power = []; Max_brake_power = []; %#ok<NASGU>
    V1_traction = []; V1_brake = []; %#ok<NASGU>
    max_accel_trac = []; max_accel_brake = []; %#ok<NASGU>
    run(rs_path);
    vars = struct('Mass', Mass, 'lambda', lambda, 'Davis', Davis, ...
        'Max_tractive_power', Max_tractive_power, 'Max_brake_power', Max_brake_power, ...
        'V1_traction', V1_traction, 'V1_brake', V1_brake, ...
        'max_accel_trac', max_accel_trac, 'max_accel_brake', max_accel_brake);
end

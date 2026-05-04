function [T_min, E_flatout] = flatout_baseline()
%FLATOUT_BASELINE  Globals-based flat-out baseline (for consistency checks).
% Builds X_flat using per-section upper bounds (same as mopso_main) then
% runs simulation_fun_CC_CR with all thresholds at max speed limit.
%
% Requires globals: vel_profile, dimension, var

    global vel_profile var dimension

    if isempty(vel_profile) || isempty(var) || isempty(dimension)
        error('flatout_baseline: globals vel_profile/var/dimension not set.');
    end

    % Build upper bound: same logic as mopso_main / nsga2_main_rl_sde
    n_sec = size(vel_profile, 1) - 1;
    ub = zeros(1, dimension);
    idx = 1;
    for s = 1:n_sec
        v_sec = vel_profile(s, 2);   % km/h speed limit for section s
        for k = 1:var(s)
            ub(idx) = v_sec;
            idx = idx + 1;
        end
    end

    X_flat = ub;
    [T_arr, E_flatout] = simulation_fun_CC_CR(X_flat);
    T_min = sum(T_arr);
end

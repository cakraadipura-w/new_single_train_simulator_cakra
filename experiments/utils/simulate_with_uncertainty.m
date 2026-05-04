function [T, E] = simulate_with_uncertainty(X, scenario)
%SIMULATE_WITH_UNCERTAINTY  Evaluate CC_CR solution under uncertain conditions.
%
% Used in Experiment E4 (Robustness Evaluation).
% Temporarily overrides global physics parameters, runs the simulation,
% then restores all globals to their original values.
%
% Inputs:
%   X        : 1×dim decision variable vector (CC_CR thresholds, km/h)
%   scenario : struct with fields:
%       .mass_factor    — multiplier for Mass (1.0=AW0, 1.27=AW2, 1.45=AW3)
%       .eta_traction   — traction efficiency (0.75 / 0.85 / 0.90)
%       .speed_reduction— fraction to reduce speed limits (0 or 0.20)
%
% Outputs:
%   T  : total running time (s)
%   E  : electrical energy (kWh) = E_mech / eta_traction

    global Mass lambda inertial_mass Davis
    global max_accel_trac max_accel_brake
    global vel_profile

    Mass_nom          = Mass;
    inertial_mass_nom = inertial_mass;
    Davis_nom         = Davis;
    max_accel_trac_nom  = max_accel_trac;
    max_accel_brake_nom = max_accel_brake;
    vel_profile_nom   = vel_profile;

    try
        mf = 1.0;
        if isfield(scenario,'mass_factor'), mf = scenario.mass_factor; end
        if abs(mf - 1.0) > 1e-9
            Mass          = Mass_nom * mf;
            inertial_mass = Mass * (1 + lambda);
            Davis(1)      = Davis_nom(1) * mf;
            max_accel_trac  = max_accel_trac_nom  / mf;
            max_accel_brake = max_accel_brake_nom / mf;
        end

        sr = 0;
        if isfield(scenario,'speed_reduction'), sr = scenario.speed_reduction; end
        if sr > 0
            n_sec   = size(vel_profile, 1);
            i_start = max(1, floor(n_sec/3));
            i_end   = min(n_sec, ceil(2*n_sec/3));
            vel_profile(i_start:i_end, 2) = ...
                vel_profile_nom(i_start:i_end, 2) * (1 - sr);
        end

        [T_arr, E_mech] = simulation_fun_CC_CR(X);
        T = sum(T_arr);

        eta = 0.85;
        if isfield(scenario,'eta_traction'), eta = scenario.eta_traction; end
        E = E_mech / eta;

    catch ME
        T = 1e6; E = 1e6;
        warning('simulate_with_uncertainty: %s', ME.message);
    end

    % Always restore globals
    Mass          = Mass_nom;
    inertial_mass = inertial_mass_nom;
    Davis         = Davis_nom;
    max_accel_trac  = max_accel_trac_nom;
    max_accel_brake = max_accel_brake_nom;
    vel_profile   = vel_profile_nom;
end

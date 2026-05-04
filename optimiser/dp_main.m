function [pop, GD, SP] = dp_main(vel_profile)
%DP_MAIN  Dynamic-Programming train energy optimiser (multi-weight Pareto sweep).
%
% Computes the exact Pareto-optimal trade-off curve between journey time and
% energy consumption using backward DP on a discretised (position, speed) grid.
%
% Method:
%   For each weight w ∈ [0,1] the scalar cost per step is:
%       c_w(a) = w * Δt  +  (1-w) * ΔE
%   Backward induction from terminal station gives the optimal velocity
%   trajectory.  Sweeping w produces the Pareto frontier.
%
% State space:
%   position  : s_k = k * DS  metres,  k = 0..K
%   speed     : v_j = j * DV  m/s,     j = 0..J
%
% Requires globals: vel_profile, gradient, dimension, Mass (for physics),
%   inertial_mass, Davis, Max_tractive_power, max_accel_trac, max_accel_brake
%   (all set by setup_project)
%
% Returns:
%   pop  : P × (dim+4) [zeros(dim) | T | E | rank | crowd]
%          X columns are zero (DP does not use CC_CR decision variables)
%   GD, SP : NaN vectors (not applicable)

    global vel_profile gradient dimension
    global inertial_mass Davis Max_tractive_power max_accel_trac max_accel_brake
    global pop_size iterations

    % ---- DP grid resolution ----
    DS = 5;        % position step (m)
    DV = 1.0;      % speed step (m/s)
    BRAKE_REGEN = 0.6;  % regenerative braking efficiency

    % ---- route discretisation ----
    S_max = max(vel_profile(:,1)) * 1000;  % m
    K     = floor(S_max / DS);
    s_grid = (0:K) * DS;                   % positions (m)

    % Speed limit at each grid position (interpolate vel_profile)
    vp_dist  = vel_profile(:,1) * 1000;    % km→m
    vp_speed = vel_profile(:,2) / 3.6;     % km/h→m/s
    vlim_grid = interp1(vp_dist, vp_speed, s_grid, 'previous', vp_speed(end));

    % Gradient at each grid position (‰ → m/s² equivalent: g_acc = g * grade/1000)
    if ~isempty(gradient) && size(gradient,1) > 1
        gr_dist  = gradient(:,1) * 1000;   % km→m
        gr_grade = gradient(:,2);          % ‰
        g_acc_grid = interp1(gr_dist, gr_grade * 9.81/1000, s_grid, 'linear', 0);
    else
        g_acc_grid = zeros(1, K+1);
    end

    V_max_ms = max(vlim_grid);
    J        = ceil(V_max_ms / DV);
    v_grid   = (0:J) * DV;                 % speed states (m/s)

    % ---- weight sweep for Pareto front ----
    w_vals = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0];
    n_w    = length(w_vals);

    T_out = zeros(n_w, 1);
    E_out = zeros(n_w, 1);

    fprintf('[DP] grid: K=%d positions, J=%d speeds | %d weight sweeps\n', K, J, n_w);

    for wi = 1:n_w
        w = w_vals(wi);

        % ---- backward DP ----
        INF_COST = 1e12;
        cost_to_go = INF_COST * ones(K+1, J+1);  % (position, speed) → cost

        % Terminal condition: must arrive at speed = 0
        cost_to_go(K+1, 1) = 0;           % j=0 → v=0 m/s (index 1)

        policy = zeros(K+1, J+1, 'int16'); % optimal next-speed index

        for k = K:-1:1
            v_lim_next = vlim_grid(k+1);   % speed limit at next position

            for j = 1:J+1
                v_cur = v_grid(j);

                % --- compute reachable next speeds ---
                % Max/min speed achievable in one DS step:
                a_trac = max_accel_trac - g_acc_grid(k) ...
                       - (Davis(1) + Davis(2)*v_cur + Davis(3)*v_cur^2) / inertial_mass;
                a_trac = max(0, min(a_trac, max_accel_trac));
                a_brk  = max_accel_brake + g_acc_grid(k) ...
                       - (Davis(1) + Davis(2)*v_cur + Davis(3)*v_cur^2) / inertial_mass;
                a_brk  = max(0, min(a_brk, max_accel_brake));

                v_max_next = min(v_lim_next, sqrt(max(0, v_cur^2 + 2*a_trac*DS)));
                v_min_next = max(0,           sqrt(max(0, v_cur^2 - 2*a_brk*DS)));

                j_min = max(1,   floor(v_min_next / DV) + 1);
                j_max = min(J+1, ceil(v_max_next  / DV) + 1);

                best_cost = INF_COST;
                best_j    = 0;

                for jn = j_min:j_max
                    v_next = v_grid(jn);
                    if v_next > v_lim_next + 1e-6, continue; end
                    if cost_to_go(k+1, jn) >= INF_COST, continue; end

                    % --- step cost ---
                    v_avg = (v_cur + v_next) / 2;
                    if v_avg < 1e-3
                        dT = sqrt(2 * DS / max(a_trac, 0.01));
                    else
                        dT = DS / v_avg;
                    end

                    a_step = (v_next^2 - v_cur^2) / (2*DS);
                    resist = (Davis(1) + Davis(2)*v_avg + Davis(3)*v_avg^2) / inertial_mass;
                    net_a  = a_step + resist - g_acc_grid(k);

                    if net_a > 0
                        % traction
                        P_trac = min(net_a * inertial_mass * 1000 * v_avg, Max_tractive_power);
                        dE = P_trac * dT / (3.6e6);  % kWh
                    else
                        % braking (regenerative)
                        F_brk = abs(net_a) * inertial_mass * 1000;
                        dE = -BRAKE_REGEN * F_brk * DS / (3.6e6);  % kWh (negative = recovered)
                        dE = max(dE, -1e-4);  % can't recover more than consumed
                    end

                    c = w * dT + (1 - w) * max(0, dE) + cost_to_go(k+1, jn);
                    if c < best_cost
                        best_cost = c;
                        best_j    = jn;
                    end
                end

                if best_j > 0
                    cost_to_go(k, j) = best_cost;
                    policy(k, j)     = best_j;
                end
            end
        end

        % ---- reconstruct trajectory from (k=1, j=1) = start at rest ----
        T_total = 0; E_total = 0;
        j_cur = 1;  % start at speed = 0
        feasible = true;
        for k = 1:K
            jn = policy(k, j_cur);
            if jn == 0, feasible = false; break; end

            v_cur  = v_grid(j_cur);
            v_next = v_grid(jn);
            v_avg  = (v_cur + v_next) / 2;

            if v_avg < 1e-3
                a_trac = max_accel_trac;
                dT = sqrt(2 * DS / max(a_trac, 0.01));
            else
                dT = DS / v_avg;
            end

            a_step = (v_next^2 - v_cur^2) / (2*DS);
            resist = (Davis(1) + Davis(2)*v_avg + Davis(3)*v_avg^2) / inertial_mass;
            net_a  = a_step + resist - g_acc_grid(k);

            if net_a > 0
                P_trac = min(net_a * inertial_mass * 1000 * v_avg, Max_tractive_power);
                dE = P_trac * dT / 3.6e6;
            else
                F_brk = abs(net_a) * inertial_mass * 1000;
                dE = -BRAKE_REGEN * F_brk * DS / 3.6e6;
                dE = max(dE, -1e-4);
            end

            T_total = T_total + dT;
            E_total = E_total + dE;
            j_cur = jn;
        end

        if feasible && j_cur == 1
            T_out(wi) = T_total;
            E_out(wi) = max(0, E_total);
        else
            T_out(wi) = NaN;
            E_out(wi) = NaN;
            fprintf('[DP] w=%.2f: infeasible trajectory\n', w);
        end

        fprintf('[DP] w=%.2f | T=%.1f s | E=%.4f kWh\n', w, T_out(wi), E_out(wi));
    end

    % ---- assemble population matrix ----
    valid = ~isnan(T_out) & ~isnan(E_out);
    T_v = T_out(valid); E_v = E_out(valid);
    n_valid = sum(valid);

    if n_valid == 0
        warning('dp_main: no feasible DP solutions found.');
        pop = zeros(0, dimension+4);
        GD = NaN; SP = NaN;
        return;
    end

    % NSGA-II style output: [X (zeros) | T | E | rank | crowd]
    F_mat = [T_v, E_v];
    rank_v  = pareto_rank(F_mat);
    crowd_v = pareto_crowd(F_mat);

    pop = [zeros(n_valid, dimension), T_v, E_v, rank_v, crowd_v];

    GD = NaN(1, iterations);
    SP = NaN(1, iterations);

    fprintf('[DP] done | %d feasible solutions on approximate Pareto front\n', n_valid);
end

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function rank = pareto_rank(F)
    n = size(F,1); rank = ones(n,1);
    for i = 1:n
        for j = 1:n
            if i~=j && all(F(j,:)<=F(i,:)) && any(F(j,:)<F(i,:))
                rank(i) = rank(i) + 1;
            end
        end
    end
end

function cd = pareto_crowd(F)
    n = size(F,1);
    if n <= 2, cd = Inf(n,1); return; end
    cd = zeros(n,1);
    for m = 1:size(F,2)
        [~,idx] = sort(F(:,m));
        cd(idx(1)) = Inf; cd(idx(end)) = Inf;
        fr = F(idx(end),m) - F(idx(1),m) + eps;
        for i = 2:n-1
            cd(idx(i)) = cd(idx(i)) + (F(idx(i+1),m)-F(idx(i-1),m))/fr;
        end
    end
end

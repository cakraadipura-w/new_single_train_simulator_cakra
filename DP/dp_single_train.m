function [E_kWh, v_opt, t_opt, s_opt] = dp_single_train(seg_file, train, T_target, varargin)
% DP_SINGLE_TRAIN  Dynamic Programming optimasi profil kecepatan satu segmen inter-stasiun.
%
% Inputs:
%   seg_file : string, path ke file .mat segmen (berisi vel_profile, gradient, station_info)
%   train    : struct dengan field:
%               Mass (t), lambda, inertial_mass (t), Davis [A B C] (kN!),
%               max_speed_kmh, V1_traction (m/s), V2_traction (m/s),
%               Max_tractive_power (W), V1_brake (m/s), V2_brake (m/s),
%               Max_brake_power (W), max_accel_trac (m/s^2), max_accel_brake (m/s^2),
%               gravity (m/s^2)
%   T_target : scalar, target running time (s)
%   Optional (Name,Value pairs):
%       'dx' : distance step (m), default 10
%       'dv' : velocity step (m/s), default 1
%       'lambda_tol' : tolerance waktu (s) pada bisection, default 0.5
%       'max_iter' : iterasi maksimum bisection, default 30
%
% Outputs:
%   E_kWh   : total traction energy (kWh)
%   v_opt   : velocity profile (m/s) pada setiap posisi (ukuran N_s x 1)
%   t_opt   : cumulative time (s) pada setiap posisi
%   s_opt   : posisi (m) dari 0 hingga panjang segmen

%% Parse optional inputs
p = inputParser;
addParameter(p, 'dx', 10);
addParameter(p, 'dv', 1);
addParameter(p, 'lambda_tol', 0.5);
addParameter(p, 'max_iter', 30);
parse(p, varargin{:});
dx = p.Results.dx;
dv = p.Results.dv;
lambda_tol = p.Results.lambda_tol;
max_iter = p.Results.max_iter;

%% Load segment data
data = load(seg_file);
vel_profile = data.vel_profile;
gradient    = data.gradient;
station_info = data.station_info;

s_start = station_info(1,1) * 1000;
s_end   = station_info(2,1) * 1000;
L = s_end - s_start;
if L <= 0
    error('Panjang segmen tidak valid.');
end

s_vec = (0:dx:L)';
N_s = length(s_vec);

v_lim = interp1(vel_profile(:,1)*1000, vel_profile(:,2)/3.6, s_start + s_vec, 'previous');
v_lim = min(v_lim, train.max_speed_kmh/3.6);
v_lim(end) = 0;

grad_permil = interp1(gradient(:,1)*1000, gradient(:,2), s_start + s_vec, 'previous');
g_force = - (grad_permil / 1000) * train.gravity / (1 + train.lambda);

%% Diskritisasi kecepatan
v_max = max(v_lim);
v_grid = (0:dv:v_max)';
N_v = length(v_grid);

%% Precompute acceleration limits
R_acc = @(v) (train.Davis(1) + train.Davis(2)*v + train.Davis(3)*v.^2) / train.inertial_mass;

a_max_trac = zeros(N_s, N_v);
a_max_brake = zeros(N_s, N_v);
for i = 1:N_s
    for j = 1:N_v
        v = v_grid(j);
        a_trac = traction_acc(v, train);
        a_brake = brake_acc(v, train);
        a_max_trac(i,j) = a_trac - R_acc(v) + g_force(i);
        a_max_brake(i,j) = a_brake + R_acc(v) - g_force(i);
    end
end

%% Fungsi DP dengan lambda
    function [E_total, T_total, v_prof] = dp_lambda(lambda)
        cost = inf(N_s, N_v);
        next_v = zeros(N_s, N_v);
        cost(N_s, 1) = 0;
        
        for i = N_s-1:-1:1
            for j = 1:N_v
                v_cur = v_grid(j);
                if v_cur > v_lim(i) + 1e-6
                    continue;
                end
                best_cost = inf;
                best_k = 0;
                for k = 1:N_v
                    v_next = v_grid(k);
                    if v_next > v_lim(i+1) + 1e-6
                        continue;
                    end
                    a_req = (v_next^2 - v_cur^2) / (2*dx);
                    if a_req > a_max_trac(i,j) + 1e-6 || a_req < -a_max_brake(i,j) - 1e-6
                        continue;
                    end
                    if v_cur + v_next < 1e-6
                        continue;
                    end
                    dt = 2*dx / (v_cur + v_next);
                    if a_req > 0
                        F_per_m = a_req + R_acc(v_cur) - g_force(i);
                        if F_per_m > 0
                            E_step = (F_per_m * train.Mass * 1000) * dx;
                        else
                            E_step = 0;
                        end
                    else
                        E_step = 0;
                    end
                    step_cost = E_step + lambda * dt;
                    total_cost = step_cost + cost(i+1, k);
                    if total_cost < best_cost
                        best_cost = total_cost;
                        best_k = k;
                    end
                end
                if best_k > 0
                    cost(i,j) = best_cost;
                    next_v(i,j) = best_k;
                end
            end
        end
        
        % Forward recovery
        v_prof = zeros(N_s,1);
        E_total = 0;
        T_total = 0;
        j = 1;
        for i = 1:N_s-1
            v_cur = v_grid(j);
            v_prof(i) = v_cur;
            k = next_v(i,j);
            if k == 0
                error('DP gagal di posisi %d.', i);
            end
            v_next = v_grid(k);
            a_req = (v_next^2 - v_cur^2) / (2*dx);
            dt = 2*dx / (v_cur + v_next);
            if a_req > 0
                F_per_m = a_req + R_acc(v_cur) - g_force(i);
                if F_per_m > 0
                    E_total = E_total + (F_per_m * train.Mass * 1000) * dx;
                end
            end
            T_total = T_total + dt;
            j = k;
        end
        v_prof(N_s) = 0;
    end

%% Bisection untuk mencari lambda yang menghasilkan T_target (dikoreksi)
[~, T0] = dp_lambda(0);

if abs(T0 - T_target) < lambda_tol
    lambda_opt = 0;
elseif T0 > T_target
    % Perlu lambda positif (penalti) agar lebih cepat
    lambda_low = 0;
    lambda_high = 1;
    while true
        [~, T_test] = dp_lambda(lambda_high);
        if T_test < T_target
            break;
        end
        lambda_high = lambda_high * 2;
        if lambda_high > 1e8
            error('Tidak dapat menemukan lambda yang cukup besar.');
        end
    end
else
    % Perlu lambda negatif (insentif) agar lebih lambat
    lambda_low = -1;
    lambda_high = 0;
    while true
        [~, T_test] = dp_lambda(lambda_low);
        if T_test > T_target
            break;
        end
        lambda_low = lambda_low * 2;
        if lambda_low < -1e8
            error('Tidak dapat menemukan lambda negatif yang cukup besar.');
        end
    end
end

if abs(T0 - T_target) >= lambda_tol
    for iter = 1:max_iter
        lambda_mid = (lambda_low + lambda_high)/2;
        [~, T_mid] = dp_lambda(lambda_mid);
        if abs(T_mid - T_target) < lambda_tol
            break;
        end
        % T menurun terhadap lambda. 
        % Jika T_mid > T_target, lambda_mid terlalu kecil → naikkan batas bawah
        % Jika T_mid < T_target, lambda_mid terlalu besar → turunkan batas atas
        if T_mid > T_target
            lambda_low = lambda_mid;
        else
            lambda_high = lambda_mid;
        end
    end
    lambda_opt = (lambda_low + lambda_high)/2;
else
    lambda_opt = 0;
end

[E_joule, T_final, v_prof] = dp_lambda(lambda_opt);

E_kWh = E_joule / (1000 * 3600);

% Output
v_opt = v_prof;
t_opt = cumsum([0; 2*dx./(v_opt(1:end-1) + v_opt(2:end))]);
t_opt = t_opt(:);
s_opt = s_vec;

end

%% ==== SUBFUNCTIONS ====
function a = traction_acc(v, train)
if v <= train.V1_traction
    F = train.Max_tractive_power / train.V1_traction;
elseif v <= train.V2_traction
    F = train.Max_tractive_power / v;
else
    F = train.Max_tractive_power * train.V2_traction / v^2;
end
a = F / (train.inertial_mass * 1000);
if a > train.max_accel_trac
    a = train.max_accel_trac;
end
end

function a = brake_acc(v, train)
if v < 1e-3
    a = train.max_accel_brake;
    return;
end
if v <= train.V1_brake
    F = train.Max_brake_power / train.V1_brake;
elseif v <= train.V2_brake
    F = train.Max_brake_power / v;
else
    F = train.Max_brake_power * train.V2_brake / v^2;
end
a = F / (train.inertial_mass * 1000);
if a > train.max_accel_brake
    a = train.max_accel_brake;
end
end
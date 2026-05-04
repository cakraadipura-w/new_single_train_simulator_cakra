%% =========================================================================
% TEST_DP_ONE_SEGMENT_FINAL_V2.m
% Perbaikan fundamental: gaya traksi/rem dihitung dengan benar (kN -> N)
% Mode: MaxTraction, Coasting, Cruising, MaxBrake.
% =========================================================================
clear; clc; close all;

%% 1. Load Rolling Stock & Set Global (seperti baseline)
fprintf('Loading rolling stock...\n');
run('rollingstock_Guangzhou_L7.m');  % Semua parameter global terdefinisi

% Tambahan: pastikan inertial mass dalam kg
global Mass lambda inertial_mass_kg
inertial_mass_kg = Mass * (1 + lambda) * 1000;  % kg

% Definisikan batas percepatan (m/s²) yang belum ada di rollingstock
a_max = max_accel_trac;      % 1.0 (dari rollingstock)
a_min = -max_accel_brake;    % -1.0 (dari rollingstock), atau bisa -1.5 jika ingin lebih longgar

%% 2. Load Data Segmen IS01
data = load('Guangzhou_Line7_IS01_0.000-1.120km.mat');
if isfield(data, 'vel_profile')
    vp_raw = data.vel_profile;
else
    vp_raw = data.speed_limit;
end
if isfield(data, 'gradient')
    grad_raw = data.gradient;
else
    grad_raw = data.slope;
end

vel_lim = [vp_raw(:,1)*1000, vp_raw(:,2)/3.6];  % [m, m/s]
grad   = [grad_raw(:,1)*1000, grad_raw(:,2)];    % [m, ‰]
total_dist = vel_lim(end,1);

Target_T = 129;  % detik (jadwal IS01)

%% 3. Setup Diskritisasi
dx = 10;            % m
dv = 0.5;           % m/s
x_grid = 0 : dx : total_dist;
v_max = max(vel_lim(:,2)) + 1;
v_grid = 0 : dv : v_max;
epsilon_v = 0.5;    % toleransi berhenti (m/s)
N_stages = length(x_grid);
Nv = length(v_grid);

% Mode kontrol (indeks 1..4):
% 1: MaxTraction (dengan batas a_max)
% 2: Coasting
% 3: Cruising (hanya jika feasible, untuk speed-holding)
% 4: MaxBrake
N_modes = 4;

fprintf('Grid: %d stages, %d speeds, %d modes\n', N_stages, Nv, N_modes);

%% 4. Fungsi Pembantu Gaya & Percepatan (DIPERBAIKI)
function [F, a] = compute_accel(v, mode, v_lim_cur, F_grade, R_davis, inertial_mass_kg, a_max, a_min, ...
    Max_tractive_power, V1_traction, V2_traction, ...
    Max_brake_power, V1_brake, V2_brake)

    % NOTE: 
    % - Max_tractive_power = 3716.8 * 1000 W
    % - V1_traction = 46.3 / 3.6 = 12.86 m/s
    % - V2_traction = 80 / 3.6 = 22.22 m/s
    % - F_trac_max pada kecepatan rendah = P / V1 = 3716800 / 12.86 = 289,000 N = 289 kN
    %   Ini cocok dengan data: max_traction = 289 kN

    switch mode
        case 1  % MaxTraction (dgn batas a_max)
            % Hitung gaya traksi maksimum yang tersedia (N)
            if v < V1_traction
                F_trac_max = Max_tractive_power / V1_traction;  % ~289 kN
            elseif v <= V2_traction
                F_trac_max = Max_tractive_power / v;
            else
                F_trac_max = Max_tractive_power * V2_traction / v^2;
            end
            % Hitung percepatan jika traksi penuh
            a_full = (F_trac_max - R_davis - F_grade) / inertial_mass_kg;
            if a_full > a_max
                % Batasi percepatan ke a_max -> hitung gaya yang diperlukan
                F_req = inertial_mass_kg * a_max + R_davis + F_grade;
                F = min(F_req, F_trac_max);  % tidak boleh melebihi kemampuan motor
            else
                F = F_trac_max;
            end
            a = (F - R_davis - F_grade) / inertial_mass_kg;
            
        case 2  % Coasting
            F = 0;
            a = (-R_davis - F_grade) / inertial_mass_kg;
            
        case 3  % Cruising (a=0)
            F_req = R_davis + F_grade;
            if F_req >= 0
                % Perlu traksi, cek apakah cukup
                if v < V1_traction
                    F_max_avail = Max_tractive_power / V1_traction;
                elseif v <= V2_traction
                    F_max_avail = Max_tractive_power / v;
                else
                    F_max_avail = Max_tractive_power * V2_traction / v^2;
                end
                if F_req <= F_max_avail
                    F = F_req;
                    a = 0;
                else
                    % Tidak mampu cruising, fallback ke coasting
                    F = 0;
                    a = (-R_davis - F_grade) / inertial_mass_kg;
                end
            else
                % Butuh pengereman untuk cruising, fallback coasting
                F = 0;
                a = (-R_davis - F_grade) / inertial_mass_kg;
            end
            
        case 4  % MaxBrake
            if v > 1
                if v <= V1_brake
                    F_brake_max = Max_brake_power / V1_brake;  % ~352 kN
                elseif v <= V2_brake
                    F_brake_max = Max_brake_power / v;
                else
                    F_brake_max = Max_brake_power * V2_brake / v^2;
                end
            else
                F_brake_max = Max_brake_power / V1_brake;
            end
            F = -F_brake_max;  % negatif karena pengereman
            a = (F - R_davis - F_grade) / inertial_mass_kg;
            % Batasi deselerasi jika perlu
            if abs(a) > abs(a_min)
                % Kurangi gaya rem
                F_req = inertial_mass_kg * a_min + R_davis + F_grade;
                F = max(F_req, -F_brake_max);
                a = (F - R_davis - F_grade) / inertial_mass_kg;
            end
    end
end

%% 5. DP Solver
function [T_total, E_total_kWh, v_opt, t_opt, x_opt] = solve_dp(lambda, ...
    vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
    inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake, ...
    Mass, gravity, Davis)

    N_stages = length(x_grid);
    Nv = length(v_grid);
    N_modes = 4;

    J = inf(Nv, N_stages);
    policy = zeros(Nv, N_stages-1, 'int8');

    % Terminal condition: semua state dengan v <= epsilon_v valid, biaya 0
    for iv = 1:Nv
        if v_grid(iv) <= epsilon_v
            J(iv, end) = 0;
        end
    end

    % Backward DP
    for k = N_stages-1 : -1 : 1
        x_cur = x_grid(k);
        v_lim_cur = interp1(vel_lim(:,1), vel_lim(:,2), x_cur, 'previous', 0);
        grad_permil = interp1(grad(:,1), grad(:,2), x_cur, 'previous', 0);
        F_grade = Mass * 1000 * gravity * sin(atan(grad_permil/1000));

        for iv = 1:Nv
            v = v_grid(iv);
            if v > v_lim_cur + 1e-6
                continue;
            end
            v_kmh = v * 3.6;
            R_davis = (Davis(1) + Davis(2)*v_kmh + Davis(3)*v_kmh^2) * Mass * gravity / 1000;

            % Iterasi untuk setiap mode
            for mode = 1:N_modes
                % Hitung percepatan dari mode ini
                [F_applied, a] = compute_accel(v, mode, v_lim_cur, F_grade, R_davis, inertial_mass_kg, a_max, a_min, ...
                    Max_tractive_power, V1_traction, V2_traction, ...
                    Max_brake_power, V1_brake, V2_brake);

                % Cek batas percepatan
                if a > a_max || a < a_min
                    continue;
                end

                v_sq = v^2 + 2*a*dx;
                if v_sq < 0
                    continue;
                end
                v_next = sqrt(v_sq);

                dt = 2*dx / (v + max(v_next, 0.01));
                dE = max(0, F_applied) * dx;  % Joule

                % Mekanisme berhenti di langkah terakhir
                extra_T = 0;
                extra_E = 0;
                if k == N_stages-1 && v_next > epsilon_v
                    % Harus berhenti dari v_next ke epsilon_v
                    % Pakai perlambatan maksimum yang aman: a_min
                    dist_stop = (epsilon_v^2 - v_next^2) / (2 * a_min);
                    time_stop = (epsilon_v - v_next) / a_min;
                    if dist_stop > 0
                        extra_T = time_stop;
                        extra_E = 0;  % pengereman tidak konsumsi energi
                    else
                        continue;  % tidak mungkin berhenti
                    end
                end

                % Biaya total
                if k == N_stages-1 && v_next <= epsilon_v
                    J_ref = 0;
                elseif k < N_stages-1
                    [~, iv_next] = min(abs(v_grid - v_next));
                    J_ref = J(iv_next, k+1);
                else
                    J_ref = 0; % sudah ditangani extra
                end
                J_cand = dE + extra_E + lambda*(dt + extra_T) + J_ref;
                if J_cand < J(iv, k)
                    J(iv, k) = J_cand;
                    policy(iv, k) = mode;
                end
            end
        end
    end

    % Forward
    [~, iv] = min(abs(v_grid)); % v=0
    v_opt = zeros(1, N_stages);
    t_opt = zeros(1, N_stages);
    E_opt = zeros(1, N_stages);
    x_opt = x_grid;
    v_opt(1) = 0;
    for k = 1:N_stages-1
        mode = policy(iv, k);
        v = v_grid(iv);
        x_cur = x_grid(k);
        grad_permil = interp1(grad(:,1), grad(:,2), x_cur, 'previous', 0);
        F_grade = Mass * 1000 * gravity * sin(atan(grad_permil/1000));
        v_kmh = v*3.6;
        R_davis = (Davis(1) + Davis(2)*v_kmh + Davis(3)*v_kmh^2) * Mass * gravity / 1000;
        [F_applied, a] = compute_accel(v, mode, interp1(vel_lim(:,1), vel_lim(:,2), x_cur, 'previous', 0), ...
            F_grade, R_davis, inertial_mass_kg, a_max, a_min, ...
            Max_tractive_power, V1_traction, V2_traction, ...
            Max_brake_power, V1_brake, V2_brake);
        v_next = sqrt(v^2 + 2*a*dx);
        [~, iv] = min(abs(v_grid - v_next));
        v_opt(k+1) = v_next;
        dt = 2*dx / (v + max(v_next,0.01));
        t_opt(k+1) = t_opt(k) + dt;
        E_opt(k+1) = E_opt(k) + max(0, F_applied)*dx;
    end
    T_total = t_opt(end);
    E_total_kWh = E_opt(end) / 1000 / 3600;
end

%% 6. TES SATU LAMBDA DULU UNTUK VERIFIKASI
lam_test = 5;
[T, E, v_opt, t_opt, x_opt] = solve_dp(lam_test, ...
    vel_lim, grad, x_grid, v_grid, dx, a_max, a_min, epsilon_v, ...
    inertial_mass_kg, ...
    Max_tractive_power, V1_traction, V2_traction, Max_brake_power, V1_brake, V2_brake, ...
    Mass, gravity, Davis);

fprintf('Test lambda=%.1f -> T=%.1f s, E=%.2f kWh\n', lam_test, T, E);

% Plot
figure;
subplot(1,2,1);
plot(x_opt/1000, v_opt*3.6, 'b-', 'LineWidth', 2); hold on;
plot(vel_lim(:,1)/1000, vel_lim(:,2)*3.6, 'r--');
xlabel('Posisi (km)'); ylabel('Kecepatan (km/h)');
title(sprintf('DP Test λ=%.1f, T=%.1fs, E=%.2fkWh', lam_test, T, E));
grid on;

subplot(1,2,2);
yyaxis left;
plot(t_opt, v_opt*3.6, 'b-'); ylabel('Kecepatan (km/h)');
yyaxis right;
plot(t_opt(2:end), diff(v_opt)/dx, 'g-'); ylabel('Akselerasi (m/s^2)');
xlabel('Waktu (s)');
title('Kecepatan & Akselerasi vs Waktu'); grid on;
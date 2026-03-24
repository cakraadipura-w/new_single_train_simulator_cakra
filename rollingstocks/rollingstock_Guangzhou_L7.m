% rollingstock_Guangzhou_L7.m
% Berdasarkan data: 'Driving strategy optimization and field test on an urban rail transit system'

% 1. Vehicle Data
Mass          = 204;             % t (AW0)
lambda        = 0.08;            % Rotary allowance
inertial_mass = Mass * (1 + lambda); % Effective mass, t

% 2. Davis Coefficients
% Formula: [A B C] dalam unit kN atau sesuai standar simulator Anda
Davis = [27, 0 * 3.6, 0.0042 * 3.6^2]; % AW0

% 3. Limits
max_speed = 80;                  % km/h
Terminal_time = 15 * 60;         % Station terminal time (seconds)
notch_num = 7;                   % Number of notches
% dwell = 30;

% 4. Traction (Dikonversi ke Watt & Newton)
Max_tractive_power = 3716.8 * 1000; % Max tractive power (W)
% V1_traction dihitung dari (3716.8/289)*3.6 / 3.6 = 12.8611 m/s
V1_traction = 46.3 / 3.6;        % m/s
max_traction = (Max_tractive_power / 1000) / V1_traction; % kN

% 5. Braking (Dikonversi ke Watt & Newton)
Max_brake_power = 3911.2 * 1000;    % Max braking power (W)
% V1_brake dihitung dari (3911.2/352)*3.6 / 3.6 = 11.1114 m/s
V1_brake = 40.0009 / 3.6;        % m/s
max_brake = (Max_brake_power / 1000) / V1_brake; % kN

% --- PERHITUNGAN TURUNAN (Agar konsisten dengan main simulator) ---

% Titik Pindah Karakteristik (V1, V2, V3)
V2_traction = 80 / 3.6;          % m/s
V3_traction = 80 / 3.6;          % m/s

V2_brake = 80 / 3.6;             % m/s
V3_brake = 80 / 3.6;             % m/s

% Parameter Fisika Tambahan
gravity = 9.81;
co_fric = 0.1;

% Akselerasi Maksimum (m/s^2)
% Menggunakan massa kendaraan (Mass) sesuai logika file original
max_accel_trac = max_traction / Mass; 
max_accel_brake = max_brake / Mass;
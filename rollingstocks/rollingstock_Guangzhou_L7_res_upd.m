% rollingstock_Guangzhou_L7_updated.m
% Berdasarkan data: 'Driving strategy optimization and field test on an urban rail transit system'

% 1. Vehicle data
Mass          = 204;             % t (AW0)
rotary        = 0.08;            % lambda
inertial_mass = Mass * (1 + rotary); % effective mass, t

% 2. Davis coefficients (Dikonversi ke kN)
% Formula asli: [27 0*3.6 0.0042*3.6^2]*Mass/1000
Davis = [27, 0, 0.0042 * 3.6^2] * Mass / 1000; 

% 3. Limits
max_speed = 80;                  % km/h

% 4. Traction (Dikonversi ke Watt & Newton)
% Efisiensi engine 0.82 sudah termasuk di sini
Max_tractive_power = 3716.8 * 1000 * 0.82; 
% V1_traction asli adalah 37.97/3.6 m/s
max_traction       = Max_tractive_power / (37.97 / 3.6); 

% 5. Braking (Dikonversi ke Watt & Newton)
Max_brake_power    = 3911.2 * 1000;       
% V1_brake asli adalah 40.0009/3.6 m/s
max_brake          = Max_brake_power / (40.0009 / 3.6); 
brake_reg_efficiency = 0.80; % Asumsi efisiensi regeneratif untuk sistem Metro L7

% --- PERHITUNGAN TURUNAN (Agar konsisten dengan main simulator) ---

% Titik Pindah Karakteristik (V1, V2, V3)
V1_traction = Max_tractive_power / max_traction; 
V2_traction = max_speed / 3.6; 
V3_traction = max_speed / 3.6;

V1_brake = Max_brake_power / max_brake;
V2_brake = max_speed / 3.6;
V3_brake = max_speed / 3.6;

% Akselerasi Maksimum (m/s^2)
% a = F / m (Newton / kg)
max_accel_trac = max_traction / (Mass * 1000);
max_accel_brake = max_brake / (Mass * 1000);

gravity = 9.81;
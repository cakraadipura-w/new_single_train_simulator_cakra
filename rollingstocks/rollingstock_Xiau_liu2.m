% rollingstock_Xiau_liu_DYNAMIC.m
% Rolling stock parameters (Table 2) - Liu et al., IET ITS 2024.

% 1. Vehicle data
Mass          = 175.6;           % t
rotary        = 0.07;            % rotary allowance (lambda)
inertial_mass = Mass * (1 + rotary); % effective mass, t

% 2. Davis coefficients (kN)
Davis = [3.856, 0.267, 0.020];

% 3. Limits
max_speed = 80;                  % km/h

% 4. Traction (SUDAH DIKALI 1000 AGAR JADI WATT & NEWTON)
Max_tractive_power = 2387 * 1000;       % Watt
max_traction       = 191.7 * 1000;      % Newton

% 5. Braking (SUDAH DIKALI 1000 AGAR JADI WATT & NEWTON)
Max_brake_power    = 2377 * 1000;       % Watt
max_brake          = 142.6 * 1000;      % Newton
brake_reg_efficiency = 0.6;

% --- TAMBAHAN PENTING (MENIRU TBILISI) ---
% Menghitung V1, V2, Acc di sini agar main.m tidak salah hitung.
 n 
% Titik Pindah Gigi (V1, V2)
% V1 = Power / Force (Watt / Newton = m/s)
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

gravity=9.8;


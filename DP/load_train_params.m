function train = load_train_params(param_script)
% LOAD_TRAIN_PARAMS  Baca file .m parameter kereta, bangun struct 'train'.
%
% Input : param_script  - path ke file .m (contoh: 'rollingstock_Guangzhou_L7.m')
% Output: train         - struct dengan semua field yang diperlukan dp_single_train

% Jalankan script untuk mengeksekusi semua variabel ke workspace fungsi ini
run(param_script);

% Cek variabel yang wajib ada
required = {'Mass','lambda','inertial_mass','Davis','max_speed',...
            'Max_tractive_power','V1_traction','V2_traction',...
            'Max_brake_power','V1_brake','V2_brake','gravity'};
for i = 1:length(required)
    if ~exist(required{i},'var')
        error('Variabel ''%s'' tidak ditemukan di %s.', required{i}, param_script);
    end
end

% Bangun struct
train.Mass               = Mass;                  % t
train.lambda             = lambda;
train.inertial_mass      = inertial_mass;         % t
train.Davis              = Davis;                 % kN (A, B, C)
train.max_speed_kmh      = max_speed;             % km/h
train.V1_traction        = V1_traction;           % m/s
train.V2_traction        = V2_traction;           % m/s
train.Max_tractive_power = Max_tractive_power;    % W
train.V1_brake           = V1_brake;              % m/s
train.V2_brake           = V2_brake;              % m/s
train.Max_brake_power    = Max_brake_power;       % W
train.gravity            = gravity;

% Kalkulasi akselerasi maksimum (dari gaya konstan di kecepatan rendah)
max_trac_force  = Max_tractive_power / 1000 / V1_traction;  % kN
max_brake_force = Max_brake_power  / 1000 / V1_brake;       % kN
train.max_accel_trac  = max_trac_force  / Mass;             % m/s^2
train.max_accel_brake = max_brake_force / Mass;             % m/s^2

end
% main_run_dp.m
% ============= MAIN SCRIPT =============
% Menjalankan DP optimal speed profile untuk segmen terpilih.
% Cukup edit bagian "KONFIGURASI" di bawah.
% Log lengkap di command window + opsional simpan ke file.
%
% Dibutuhkan di path yang sama:
%   - dp_single_train.m   (template DP yang sudah diberikan sebelumnya)
%   - file .mat rute setiap segmen
%   - file rolling stock (.m) yang akan dipanggil (misal rollingstock_Guangzhou_L7.m)

clear; close all; clc;

% ==========================================
%          KONFIGURASI (EDIT SESUAI KEBUTUHAN)
% ==========================================

% 1. Pilih file rolling stock
rolling_stock_file = 'rollingstock_Guangzhou_L7.m';

% 2. Daftar segmen dan jadwal waktu tempuh (T_sched)
ALL_SEGMENTS = struct( ...
    'name',    {'IS01','IS02','IS03','IS04','IS05','IS06','IS07','IS08'}, ...
    'file',    { ...
        'Guangzhou_Line7_IS01_0.000-1.120km.mat', ...
        'Guangzhou_Line7_IS02_1.120-3.028km.mat', ...
        'Guangzhou_Line7_IS03_3.028-5.200km.mat', ...
        'Guangzhou_Line7_IS04_5.200-6.842km.mat', ...
        'Guangzhou_Line7_IS05_6.842-8.958km.mat', ...
        'Guangzhou_Line7_IS06_8.958-11.323km.mat', ...
        'Guangzhou_Line7_IS07_11.323-13.729km.mat', ...
        'Guangzhou_Line7_IS08_13.729-17.507km.mat'}, ...
    'T_sched', {129, 170, 185, 180, 185, 220, 210, 330});

% 3. Opsi DP
dp_options = struct('dx', 10, 'dv', 1, 'lambda_tol', 0.5, 'max_iter', 30);

% 4. Nama file log (string kosong -> tidak simpan log ke file)
log_file = 'dp_log.txt';

% 5. Simpan hasil ke file .mat?
save_results = false;
results_file = 'dp_results.mat';

% 6. Pilih segmen yang akan dijalankan (kosongkan [] untuk semua)
%    Bisa pakai indeks, contoh: [1, 3, 5]
%    Bisa pakai nama,   contoh: {'IS01', 'IS03', 'IS05'}
%selected_segments = [];   % [] artinya semua segmen
 selected_segments = [1];         % hanya IS01
% selected_segments = {'IS01','IS03'};   % hanya IS01 dan IS03

% ==========================================
%          PROSES UTAMA
% ==========================================

% Setup logging
if ~isempty(log_file)
    diary(log_file);
    diary on;
end

fprintf('========================================================\n');
fprintf('   DP OPTIMAL SPEED PROFILE - MAIN RUN\n');
fprintf('   Rolling stock : %s\n', rolling_stock_file);
fprintf('   Tanggal       : %s\n', datestr(now));
fprintf('========================================================\n\n');

% Load parameter kereta
fprintf('[1/4] Memuat parameter kereta... ');
try
    train = load_train_params(rolling_stock_file);
    fprintf('Berhasil.\n');
    fprintf('      Massa (AW0)        : %.1f t\n', train.Mass);
    fprintf('      Inertial mass      : %.2f t\n', train.inertial_mass);
    fprintf('      Max speed          : %.0f km/h\n', train.max_speed_kmh);
    fprintf('      Max tractive power : %.1f kW\n', train.Max_tractive_power/1000);
    fprintf('      Max brake power    : %.1f kW\n', train.Max_brake_power/1000);
    fprintf('      Max accel trac     : %.3f m/s^2\n', train.max_accel_trac);
    fprintf('      Max accel brake    : %.3f m/s^2\n', train.max_accel_brake);
catch ME
    fprintf('GAGAL.\n');
    error('Error loading train params: %s', ME.message);
end

n_seg = length(ALL_SEGMENTS);

% Filter segmen yang akan dijalankan
if ~isempty(selected_segments)
    if isnumeric(selected_segments)
        run_idx = selected_segments;
    elseif iscell(selected_segments)
        run_idx = [];
        all_names = {ALL_SEGMENTS.name};
        for j = 1:length(selected_segments)
            idx = find(strcmp(all_names, selected_segments{j}));
            if ~isempty(idx)
                run_idx(end+1) = idx;
            else
                warning('Segmen ''%s'' tidak ditemukan, dilewati.', selected_segments{j});
            end
        end
    else
        error('Format selected_segments tidak valid.');
    end
    run_idx = run_idx(run_idx >= 1 & run_idx <= n_seg);
else
    run_idx = 1:n_seg;
end

Energy_all = zeros(n_seg,1);
Time_all   = zeros(n_seg,1);
v_profiles = cell(n_seg,1);
s_profiles = cell(n_seg,1);

fprintf('\n[2/4] Menjalankan DP untuk %d segmen...\n', length(run_idx));

for i = run_idx
    seg_name = ALL_SEGMENTS(i).name;
    seg_file = ALL_SEGMENTS(i).file;
    T_tgt    = ALL_SEGMENTS(i).T_sched;
    
    fprintf('\n--- %s (%s) ---\n', seg_name, seg_file);
    fprintf('    Jadwal T_target = %.0f s\n', T_tgt);
    
    if ~exist(seg_file, 'file')
        fprintf('    ERROR: File segmen tidak ditemukan! Dilewati.\n');
        continue;
    end
    
    t_start = tic;
    try
        [E, v, t, s] = dp_single_train(seg_file, train, T_tgt, ...
            'dx', dp_options.dx, ...
            'dv', dp_options.dv, ...
            'lambda_tol', dp_options.lambda_tol, ...
            'max_iter', dp_options.max_iter);
        elapsed = toc(t_start);
        
        Energy_all(i) = E;
        Time_all(i)   = t(end);
        v_profiles{i} = v;
        s_profiles{i} = s;
        
        fprintf('    Waktu tempuh aktual  : %.2f s\n', t(end));
        fprintf('    Energi traksi        : %.4f kWh\n', E);
        fprintf('    Kecepatan maksimum   : %.1f km/h\n', max(v)*3.6);
        fprintf('    Kecepatan rata-rata  : %.1f km/h\n', (max(s)/1000) / (t(end)/3600));
        fprintf('    Waktu komputasi      : %.2f s\n', elapsed);
    catch ME
        elapsed = toc(t_start);
        fprintf('    ERROR: %s\n', ME.message);
        fprintf('    Waktu sebelum error  : %.2f s\n', elapsed);
    end
end

fprintf('\n[3/4] Rekapitulasi Hasil\n');
fprintf('----------------------------------------------\n');
fprintf('Segmen   T_target(s)   T_act(s)   Energi(kWh)\n');
fprintf('----------------------------------------------\n');
for i = run_idx
    fprintf('%-8s %9.0f   %9.2f   %11.4f\n', ...
        ALL_SEGMENTS(i).name, ...
        ALL_SEGMENTS(i).T_sched, ...
        Time_all(i), ...
        Energy_all(i));
end
fprintf('----------------------------------------------\n');
fprintf('TOTAL (terpilih)     %9.2f   %11.4f\n', sum(Time_all(run_idx)), sum(Energy_all(run_idx)));

% Simpan hasil
if save_results
    fprintf('\n[4/4] Menyimpan hasil ke %s... ', results_file);
    save(results_file, 'ALL_SEGMENTS', 'Energy_all', 'Time_all', ...
        'v_profiles', 's_profiles', 'train', 'dp_options', 'run_idx');
    fprintf('Selesai.\n');
end

% Plot profil kecepatan segmen terpilih
figure('Name', 'Optimal Speed Profiles', 'NumberTitle', 'off');
hold on;
colors = lines(length(run_idx));
leg_str = cell(length(run_idx),1);
for idx = 1:length(run_idx)
    i = run_idx(idx);
    if ~isempty(v_profiles{i})
        plot(s_profiles{i}/1000, v_profiles{i}*3.6, 'Color', colors(idx,:), 'LineWidth', 1.5);
        leg_str{idx} = sprintf('%s (%.1f kWh)', ALL_SEGMENTS(i).name, Energy_all(i));
    end
end
hold off;
grid on;
xlabel('Jarak (km)');
ylabel('Kecepatan (km/h)');
title(sprintf('Profil Kecepatan Optimal - Total Energi %.3f kWh', sum(Energy_all(run_idx))));
legend(leg_str, 'Location', 'bestoutside');
set(gcf, 'Position', [100 100 1200 600]);

fprintf('\n========================================================\n');
fprintf('   SELESAI.\n');
fprintf('========================================================\n');

% Matikan diary
if ~isempty(log_file)
    diary off;
    fprintf('Log disimpan di %s\n', log_file);
end

%% ========== FUNGSI LOKAL ==========
function train = load_train_params(param_script)
% LOAD_TRAIN_PARAMS  Baca file .m parameter kereta, bangun struct 'train'.
    run(param_script);
    required = {'Mass','lambda','inertial_mass','Davis','max_speed',...
                'Max_tractive_power','V1_traction','V2_traction',...
                'Max_brake_power','V1_brake','V2_brake','gravity'};
    for i = 1:length(required)
        if ~exist(required{i},'var')
            error('Variabel ''%s'' tidak ditemukan di %s.', required{i}, param_script);
        end
    end
    train.Mass               = Mass;
    train.lambda             = lambda;
    train.inertial_mass      = inertial_mass;
    train.Davis              = Davis;
    train.max_speed_kmh      = max_speed;
    train.V1_traction        = V1_traction;
    train.V2_traction        = V2_traction;
    train.Max_tractive_power = Max_tractive_power;
    train.V1_brake           = V1_brake;
    train.V2_brake           = V2_brake;
    train.Max_brake_power    = Max_brake_power;
    train.gravity            = gravity;
    max_trac_force  = Max_tractive_power / 1000 / V1_traction;
    max_brake_force = Max_brake_power  / 1000 / V1_brake;
    train.max_accel_trac  = max_trac_force  / Mass;
    train.max_accel_brake = max_brake_force / Mass;
end
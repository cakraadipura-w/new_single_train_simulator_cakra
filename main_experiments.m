%% main_experiments.m — Master Orchestrator (E1–E4, configurable segments)
%
% =========================================================================
%  CARA PAKAI — edit bagian USER CONFIG di bawah, lalu jalankan file ini.
% =========================================================================
%
%  1. Pilih intersection yang mau dijalankan  (SEGMENTS_TO_RUN)
%  2. Pilih eksperimen mana yang aktif        (RUN_E1..E4)
%  3. Jalankan:  main_experiments
%
%  Hasil disimpan ke ./experiment_results/<IS_name>/

clc; clear; clear global;
addpath(genpath(pwd));

%% =========================================================================
%% ===== USER CONFIG — EDIT DI SINI ========================================
%% =========================================================================

% --- Pilih intersection yang mau dijalankan ---
% Opsi:
%   'all'               → semua 8 intersection (IS01..IS08)
%   {'IS04'}            → satu intersection
%   {'IS04','IS08'}     → beberapa intersection
%   {'IS01','IS02','IS03','IS04','IS05','IS06','IS07','IS08'}  → semua manual
SEGMENTS_TO_RUN = {'IS07'};

% --- Pilih eksperimen yang mau dijalankan ---
RUN_BASELINE = true;    % hitung T_min dan E_baseline tiap segment
RUN_E1       = true;    % Ablation Study       (4 config × 30 runs)
RUN_E2       = false;   % Tight Time Sweep     (8 slack levels × 30 runs)
RUN_E3       = false;   % Strong Benchmarking  (6 methods × 30 runs)
RUN_E4       = false;   % Robustness           (30 runs × 18 scenarios)
RUN_E5       = false;   % LLM Advisor Study    (A/B/C/D/C_LLM/D_LLM)

% --- Jumlah runs (turunkan untuk quick test) ---
N_RUNS = 1;   % 5 = quick test, 10 = medium, 30 = full experiment

%% =========================================================================
%% ===== DEFINISI SEMUA 8 SEGMENT (jangan diubah) ==========================
%% =========================================================================
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
    'T_sched', {130, 170, 185, 180, 185, 220, 210, 330} );

RS_FILE = 'rollingstock_Guangzhou_L7.m';
rs_path = which(RS_FILE);
assert(~isempty(rs_path), 'Rolling stock file tidak ditemukan: %s', RS_FILE);

% --- Resolve daftar segment yang dipilih ---
if ~iscell(SEGMENTS_TO_RUN)
    SEGMENTS_TO_RUN = {ALL_SEGMENTS.name};
end
% Validasi
for k = 1:numel(SEGMENTS_TO_RUN)
    assert(any(strcmp({ALL_SEGMENTS.name}, SEGMENTS_TO_RUN{k})), ...
        'Segment tidak dikenal: "%s". Gunakan IS01..IS08.', SEGMENTS_TO_RUN{k});
end

fprintf('============================================================\n');
fprintf('  main_experiments.m\n');
fprintf('  Segments : %s\n', strjoin(SEGMENTS_TO_RUN, ', '));
fprintf('  N_RUNS   : %d\n', N_RUNS);
fprintf('  E1=%d  E2=%d  E3=%d  E4=%d  E5=%d\n', RUN_E1,RUN_E2,RUN_E3,RUN_E4,RUN_E5);
fprintf('============================================================\n\n');

%% =========================================================================
%% ===== LOOP ATAS SETIAP SEGMENT YANG DIPILIH =============================
%% =========================================================================
for seg_k = 1:numel(SEGMENTS_TO_RUN)

    seg_name = SEGMENTS_TO_RUN{seg_k};
    seg_idx  = find(strcmp({ALL_SEGMENTS.name}, seg_name));
    ACTIVE_SEG = ALL_SEGMENTS(seg_idx);          % dipakai run_E* scripts
    ACTIVE_RS  = RS_FILE;
    ACTIVE_NRUNS = N_RUNS;

    % Output folder per-segment
    OUT_DIR = fullfile(pwd, 'experiment_results', seg_name);
    if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

    fprintf('\n####################################################\n');
    fprintf('  Segment: %s | T_sched=%d s\n', seg_name, ACTIVE_SEG.T_sched);
    fprintf('  Output : %s\n', OUT_DIR);
    fprintf('####################################################\n\n');

    %% ----- BASELINE -----
    T_min_seg   = NaN;
    E_base_seg  = NaN;
    if RUN_BASELINE
        fprintf('--- BASELINE: %s ---\n', seg_name);
        t0_bl = tic;
        rp = which(ACTIVE_SEG.file);
        if isempty(rp)
            fprintf('[SKIP] Route file tidak ditemukan: %s\n', ACTIVE_SEG.file);
        else
            try
                [T_min_seg, E_base_seg] = flat_out_baseline_rk4(rp, rs_path);
                save(fullfile(OUT_DIR,'baseline.mat'), 'T_min_seg','E_base_seg','ACTIVE_SEG');
                fprintf('T_min=%.2f s | E_base=%.4f kWh | %.1f s\n\n', ...
                    T_min_seg, E_base_seg, toc(t0_bl));
            catch ME
                fprintf('[ERROR] Baseline %s: %s\n', seg_name, ME.message);
            end
        end
    end

    %% ----- E1: ABLATION STUDY -----
    if RUN_E1
        fprintf('--- E1: Ablation Study | %s ---\n', seg_name);
        t0_E1 = tic;
        try
            run_E1_ablation;
        catch ME
            fprintf('[main] E1 ERROR (%s): %s\n', seg_name, ME.message);
            save(fullfile(OUT_DIR,'E1_error.mat'),'ME');
        end
        fprintf('[E1/%s] %.1f min\n\n', seg_name, toc(t0_E1)/60);
        keep_vars;
    end

    %% ----- E2: TIGHT TIME SWEEP -----
    if RUN_E2
        fprintf('--- E2: Tight Time Sweep | %s ---\n', seg_name); %#ok<UNRCH>
        t0_E2 = tic;
        try
            run_E2_tight_sweep;
        catch ME
            fprintf('[main] E2 ERROR (%s): %s\n', seg_name, ME.message);
            save(fullfile(OUT_DIR,'E2_error.mat'),'ME');
        end
        fprintf('[E2/%s] %.1f min\n\n', seg_name, toc(t0_E2)/60);
        keep_vars;
    end

    %% ----- E3: STRONG BENCHMARKING -----
    if RUN_E3
        fprintf('--- E3: Benchmarking | %s ---\n', seg_name); %#ok<UNRCH>
        t0_E3 = tic;
        try
            run_E3_benchmarking;
        catch ME
            fprintf('[main] E3 ERROR (%s): %s\n', seg_name, ME.message);
            save(fullfile(OUT_DIR,'E3_error.mat'),'ME');
        end
        fprintf('[E3/%s] %.1f min\n\n', seg_name, toc(t0_E3)/60);
        keep_vars;
    end

    %% ----- E4: ROBUSTNESS -----
    if RUN_E4
        fprintf('--- E4: Robustness | %s ---\n', seg_name); %#ok<UNRCH>
        t0_E4 = tic;
        try
            run_E4_robustness;
        catch ME
            fprintf('[main] E4 ERROR (%s): %s\n', seg_name, ME.message);
            save(fullfile(OUT_DIR,'E4_error.mat'),'ME');
        end
        fprintf('[E4/%s] %.1f min\n\n', seg_name, toc(t0_E4)/60);
        keep_vars;
    end

    %% ----- E5: LLM ADVISOR STUDY -----
    if RUN_E5
        fprintf('--- E5: LLM Advisor Study | %s ---\n', seg_name);
        t0_E5 = tic;
        try
            run_E5_llm_advisor;
        catch ME
            fprintf('[main] E5 ERROR (%s): %s\n', seg_name, ME.message);
            save(fullfile(OUT_DIR,'E5_error.mat'),'ME');
        end
        fprintf('[E5/%s] %.1f min\n\n', seg_name, toc(t0_E5)/60);
        keep_vars;
    end

end  % end segment loop

%% =========================================================================
%% ===== SUMMARY AKHIR ======================================================
%% =========================================================================
fprintf('\n============================================================\n');
fprintf('  SEMUA SELESAI\n');
for seg_k = 1:numel(SEGMENTS_TO_RUN)
    seg_name = SEGMENTS_TO_RUN{seg_k};
    d = fullfile(pwd,'experiment_results',seg_name);
    csv_n = numel(dir(fullfile(d,'*.csv')));
    png_n = numel(dir(fullfile(d,'*.png')));
    fprintf('  %s → %d CSV, %d PNG  (%s)\n', seg_name, csv_n, png_n, d);
end
fprintf('============================================================\n');

%% =========================================================================
%% LOCAL HELPER
function keep_vars
% Bersihkan workspace tapi pertahankan variabel kontrol orchestrator
    evalin('base', ['clearvars -except ' ...
        'ALL_SEGMENTS SEGMENTS_TO_RUN RS_FILE rs_path ' ...
        'RUN_BASELINE RUN_E1 RUN_E2 RUN_E3 RUN_E4 RUN_E5 N_RUNS ' ...
        'seg_k seg_name seg_idx ACTIVE_SEG ACTIVE_RS ACTIVE_NRUNS ' ...
        'OUT_DIR T_min_seg E_base_seg']);
end

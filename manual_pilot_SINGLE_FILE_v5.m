clc; clear; close all;

%% ==========================================================
%% 1. INPUT MANUAL ANDA DI SINI
%% ==========================================================
% Masukkan angka kecepatan untuk SETIAP section.
% Script ini akan otomatis mengecek apakah Section itu butuh 1, 2, atau 3 angka.
% (Pastikan jumlah baris input sama dengan jumlah section di rute!)

USER_INPUT_VALUES = { ...
    [35], ...         % Section 1 (Biasanya Mode CC -> 2 angka) here is the more important data
    [49, 47, 58], ...     % Section 2 (Biasanya Mode CR -> 3 angka)
    [47, 35], ...         % Section 3 (Biasanya Mode CC -> 2 angka)
    [42, 20] ...          % Section 4 (Biasanya Mode CC -> 2 angka)
};

%% ==========================================================
%% 2. LOAD DATA & GLOBAL VARS
%% ==========================================================
global pop_size iterations dimension;
global vel_profile station_info gradient terminal;
global Mass lambda inertial_mass Davis;
global V1_traction V2_traction V3_traction V1_brake V2_brake V3_brake;
global Max_tractive_power Max_brake_power max_speed co_fric gravity;
global max_traction max_brake max_accel_trac max_accel_brake;

% Global var untuk simulation_fun
global var; 

% Load Rute

%if exist('Beijing_Yizhuangline_up.mat', 'file'), load('Beijing_Yizhuangline_up.mat'); else error('File route hilang!'); end
if exist('xiau_liu_route.mat', 'file'), load('xiau_liu_route.mat'); else error('File route hilang!'); end

% Load Rolling Stock
if exist('rollingstock_Xiau_liu2.m', 'file')
    rollingstock_Xiau_liu2; 
    lambda = rotary; %gravity = 9.8;         
    if max_traction < 1000, max_traction = max_traction * 1000; end
    if max_brake < 1000, max_brake = max_brake * 1000; end
    %if Davis(1) < 100, Davis = Davis * 1000; end 
    real_mass_kg = Mass * 1000;
    
    max_accel_trac = max_traction / real_mass_kg; 
    
    max_accel_brake = max_brake / real_mass_kg;
    
    if ~exist('V1_traction', 'var')
        V1_traction = Max_tractive_power / max_traction; 
        V2_traction = max_speed / 3.6; 
        V1_brake = Max_brake_power / max_brake; V2_brake = V2_traction;
    end

% Ensure gravity is defined
if isempty(gravity) || ~isfinite(gravity)
    gravity = 9.81;
end

end

%% ==========================================================
%% ==========================================================
%% 3. STRATEGY DETECTION (IMPROVED)
%%    - segmentasi ikut breakpoint gradient
%%    - downhill -> var=3 (extra coasting DOF)
%%    - auto-normalize unit gradient (m vs km, slope vs ‰)
%% ==========================================================
fprintf('Automatic strategy detection (improved)...\n');

% Keep a raw copy of vel_profile (so we can refine repeatedly)
global vel_profile_raw
if isempty(vel_profile_raw) || ~isnumeric(vel_profile_raw) || size(vel_profile_raw,2) < 2
    vel_profile_raw = vel_profile;
end
vp_raw = vel_profile_raw;
% --- Normalize gradient units ---
% Assumptions target:
%   vp_raw(:,1) = km, vp_raw(:,2) = km/h
%   gradient(:,1) = km, gradient(:,2) = per-mille (‰)
route_km = vp_raw(end,1);

% If gradient positions look like meters -> km
if max(gradient(:,1)) > 10*max(1e-6, route_km)
    gradient(:,1) = gradient(:,1)/1000;
end
% If gradient values look like decimal slope -> per-mille
if max(abs(gradient(:,2))) < 0.2
    gradient(:,2) = gradient(:,2)*1000;
end

% --- Refine vel_profile by inserting gradient breakpoints ---
bp = unique([vp_raw(:,1); gradient(:,1)]);
bp = bp(bp>=vp_raw(1,1) & bp<=vp_raw(end,1));
bp = unique([vp_raw(1,1); bp; vp_raw(end,1)]);
bp = sort(bp);

% Inherit speed-limit from original vel_profile segment
spd = zeros(size(bp));
pos = 1;
for ii=1:numel(bp)
    while (pos < size(vp_raw,1)) && (bp(ii) >= vp_raw(pos+1,1) - 1e-12)
        pos = pos + 1;
    end
    spd(ii) = vp_raw(pos,2);
end
vel_profile = [bp(:), spd(:)];

% --- Downhill thresholds (use Davis resistance @ v_ref) ---
v_ref = 30/3.6; % m/s
Rref = Davis(1) + Davis(2)*v_ref + Davis(3)*v_ref^2; % N
m_kg = Mass*1000; % Mass given in tons in your project
min_steepdownhill = -(Rref)/(m_kg*gravity)*1000;  % ‰
min_milddownhill  = 0.5*min_steepdownhill;

% --- Build gradient segments and mark downhill sections ---
g_bp  = gradient(:,1);
g_val = gradient(:,2);
g_start = g_bp(1:end-1);
g_end   = g_bp(2:end);
g_seg   = g_val(1:end-1);

isDown = (g_seg <= min_milddownhill);

num_segments_real = size(vel_profile,1)-1;
var = 2*ones(1,num_segments_real); % default normal=2

for i=1:num_segments_real
    a = vel_profile(i,1);
    b = vel_profile(i+1,1);
    overlap = isDown & (g_start < b) & (g_end > a);
    if any(overlap)
        var(i) = 3; % downhill -> allow extra coasting threshold
    end
end

fprintf('  segments raw=%d | refined=%d | var==3 (downhill)=%d\n', ...
    size(vp_raw,1)-1, num_segments_real, nnz(var==3));

%% ==========================================================
%% 4. VALIDASI & MAPPING INPUT USER (ROBUST)
%%    NEW rule:
%%      - var==2 needs [HIGH, LOW]   (2 angka)
%%      - var==3 needs [HIGH, LOW1, LOW2] (3 angka)
%%    Tapi biar gampang, input boleh 1/2/3 angka -> script auto-complete.
%% ==========================================================
X_input = [];

fprintf('Manual Input Mapping...\n');

for i = 1:num_segments_real
    mode_req = var(i);
    vlim_i = vel_profile(i,2); % km/h
    % get user values if provided, else empty
    if i <= length(USER_INPUT_VALUES)
        vals = USER_INPUT_VALUES{i};
    else
        vals = [];
    end
    vals = vals(:)'; % row
    nvals = length(vals);

    if mode_req == 2
        % Need 2: [HIGH, LOW]  (km/h)
        if nvals == 0
            HIGH = 0.95*vlim_i; LOW = 0.70*HIGH;
            vals = [HIGH, LOW];
        elseif nvals == 1
            HIGH = vals(1); LOW = 0.70*HIGH;
            vals = [HIGH, LOW];
        else
            vals = vals(1:2);
            % enforce HIGH>=LOW
            vals = sort(vals,'descend');
        end

        % clamp to sensible range
        HIGH = min(max(vals(1), 0.2*vlim_i), vlim_i);
        LOW  = min(max(vals(2), 0.05*vlim_i), HIGH);
        vals = [HIGH, LOW];

    else % mode_req == 3
        % Need 3: [HIGH, LOW1, MID]  (km/h)
        % MID is an EARLY-coasting trigger between LOW1 and HIGH.
        if nvals == 0
            HIGH = 0.90*vlim_i;
            LOW1 = 0.60*HIGH;
            MID  = LOW1 + 0.60*(HIGH-LOW1); % between LOW1 and HIGH
            vals = [HIGH, LOW1, MID];
        elseif nvals == 1
            HIGH = vals(1);
            LOW1 = 0.60*HIGH;
            MID  = LOW1 + 0.60*(HIGH-LOW1);
            vals = [HIGH, LOW1, MID];
        elseif nvals == 2
            tmp = sort(vals(1:2),'descend');
            HIGH = tmp(1); LOW1 = tmp(2);
            if (HIGH-LOW1) < 0.12*vlim_i
                LOW1 = 0.60*HIGH;
            end
            MID  = LOW1 + 0.60*(HIGH-LOW1);
            vals = [HIGH, LOW1, MID];
        else
            tmp = sort(vals(1:2),'descend');
            HIGH = tmp(1); LOW1 = tmp(2);
            MID  = vals(3);
            if (HIGH-LOW1) < 0.12*vlim_i
                LOW1 = 0.60*HIGH;
            end
            % Clamp MID into [LOW1, HIGH]
            MID = min(max(MID, LOW1), HIGH);
            % If MID too close to HIGH, pull it down a bit (so it actually triggers earlier)
            if (HIGH-MID) < 0.08*vlim_i
                MID = LOW1 + 0.60*(HIGH-LOW1);
            end
            vals = [HIGH, LOW1, MID];
        end

        % clamp to sensible range and enforce LOW1 <= MID <= HIGH <= vlim
        HIGH = min(max(vals(1), 0.2*vlim_i), vlim_i);
        LOW1 = min(max(vals(2), 0.05*vlim_i), HIGH);
        MID  = min(max(vals(3), LOW1), HIGH);
        vals = [HIGH, LOW1, MID];
    end

    fprintf('  Section %d (var=%d, vlim=%.1f): [%s]  (HIGH-LOW1=%.1f | HIGH-MID=%.1f)\n', ...
        i, mode_req, vlim_i, num2str(vals), vals(1)-vals(2), vals(1)-vals(end));
    X_input = [X_input, vals];
end

fprintf('Total decision variables: %d\n', length(X_input));

%% 5. EKSEKUSI SIMULASI
%% ==========================================================

run = true;

if run
    tStart = tic;
    [running_inter, Total_E, s_out, v_out, vlim_out] = simulation_fun_CC_CR(X_input);
    %disp(running_inter);
    tElapsed = toc(tStart);
    
    total_time = sum(running_inter);
    
    fprintf('---------------------------------------------\n');
    fprintf('Simulation Results:\n');
    fprintf('Total of Time  : %.2f detik\n', total_time);
    fprintf('Total of Energy : %.4f kWh\n', Total_E);
    fprintf('Elapsed runtime (1 eval): %.4f detik\n', tElapsed);


%% ===== QUICK A/B CHECK (same X, but force all sections to var=2) =====
try
    var_saved = var;
    vp_saved = vel_profile;
    % force var=2, rebuild a compatible X by dropping the extra 3rd param in var=3 sections
    var = 2*ones(1, size(vel_profile,1)-1);
    X2 = [];
    k = 1;
    for ii=1:numel(var_saved)
        if var_saved(ii)==2
            X2 = [X2, X_input(k), X_input(k+1)];
            k = k + 2;
        else
            % take HIGH, LOW1 from var=3 as (cruising, coasting) proxy
            X2 = [X2, X_input(k), X_input(k+1)];
            k = k + 3;
        end
    end
    [t2,E2] = simulation_fun_CC_CR(X2);
    fprintf('\nA/B check: Forced var=2 -> Time=%.2f s | Energy=%.4f kWh\n', sum(t2), E2);
    % restore
    var = var_saved;
    vel_profile = vp_saved;
catch ME
    fprintf('\nA/B check skipped: %s\n', ME.message);
end
    fprintf('---------------------------------------------\n');
    
    % Plotting
    figure('Name', 'Single File Pilot', 'Color', 'w');
    
    subplot(2,1,1);
    plot(s_out/1000, vlim_out*3.6, '--r', 'LineWidth', 1.5); hold on;
    plot(s_out/1000, v_out*3.6, '-b', 'LineWidth', 2);
    
    % Garis batas section
    for k = 1:size(vel_profile, 1)-1
         xline(vel_profile(k,1), ':k');
    end
    
    ylabel('Speed (km/h)'); title(sprintf('Time: %.1fs, Energy: %.2f kWh', total_time, Total_E)); 
    legend('Limit', 'Actual'); grid on;
    
    subplot(2,1,2);
    grad_interp = interp1(gradient(:,1)*1000, gradient(:,2), s_out, 'previous');
    area(s_out/1000, cumtrapz(s_out, grad_interp)/1000, 'FaceColor', [0.8 0.8 0.8]);
    ylabel('Elevation (m)'); xlabel('Distance (km)'); grid on;

end
    
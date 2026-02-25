function decision_var_NO_improved()
% decision_var_NO_improved_lite.m
% Goal: keep the *benefit* of improved downhill detection, but avoid
% dimension explosion by only splitting vel_profile at "downhill state"
% transitions (instead of *every* gradient breakpoint).
%
% Output:
%   global var : values in {2,3}
%   global vel_profile : refined ONLY at downhill transitions
%
% Assumptions after normalization:
%   vel_profile(:,1) in km, vel_profile(:,2) in km/h
%   gradient(:,1) in km, gradient(:,2) in ‰

global var vel_profile gradient Davis Mass gravity
global vel_profile_raw gradient_raw

if isempty(gravity), gravity = 9.81; end

% --- basic checks ---
if isempty(vel_profile) || ~isnumeric(vel_profile) || size(vel_profile,2) < 2 || size(vel_profile,1) < 2
    var = [];
    return;
end
if isempty(gradient) || ~isnumeric(gradient) || size(gradient,2) < 2 || size(gradient,1) < 2
    var = 2*ones(1, size(vel_profile,1)-1);
    return;
end

% --- store raw once per route ---
if isempty(vel_profile_raw) || ~isnumeric(vel_profile_raw) || size(vel_profile_raw,2) < 2
    vel_profile_raw = vel_profile;
end
if isempty(gradient_raw) || ~isnumeric(gradient_raw) || size(gradient_raw,2) < 2
    gradient_raw = gradient;
end

vp_raw = vel_profile_raw;
gradient = gradient_raw;

% --- normalize units for gradient ---
route_km = vp_raw(end,1);

if max(gradient(:,1)) > 10*max(1e-6, route_km)   % probably meters
    gradient(:,1) = gradient(:,1) / 1000;
end
if max(abs(gradient(:,2))) < 0.2                % probably slope (e.g., 0.012)
    gradient(:,2) = gradient(:,2) * 1000;       % -> ‰
end

% --- downhill thresholds (same spirit as original) ---
v_ref = 30/3.6; % m/s
min_steepdownhill = -(Davis(1) + Davis(2)*v_ref + Davis(3)*v_ref^2) / Mass / gravity * 1000; % ‰
min_milddownhill  = 0.5 * min_steepdownhill; % negative

% --- gradient segments + downhill state ---
g_bp  = gradient(:,1);
g_val = gradient(:,2);
g_start = g_bp(1:end-1);
g_end   = g_bp(2:end);
g_seg   = g_val(1:end-1);

isDown = (g_seg <= min_milddownhill);

% --- pick only breakpoints where downhill state changes ---
chg = [true; (isDown(2:end) ~= isDown(1:end-1)); true];
bp_down = unique([g_start(chg(1:end-1)); g_end(chg(2:end))]);

% --- refine vel_profile only with these important points ---
vp_min = vp_raw(1,1);
vp_max = vp_raw(end,1);

bp = unique([vp_raw(:,1); bp_down(:)]);
bp = bp(bp >= vp_min & bp <= vp_max);
bp = unique([vp_min; bp; vp_max]);
bp = sort(bp);

% inherit speed limit from raw vel_profile segment
spd = zeros(size(bp));
pos = 1;
for i=1:numel(bp)
    while (pos < size(vp_raw,1)) && (bp(i) >= vp_raw(pos+1,1) - 1e-12)
        pos = pos + 1;
    end
    spd(i) = vp_raw(pos,2);
end
vel_profile = [bp(:), spd(:)];

% --- assign var per refined section ---
nSec = size(vel_profile,1) - 1;
var  = 2 * ones(1, nSec);

for i=1:nSec
    a = vel_profile(i,1);
    b = vel_profile(i+1,1);
    overlap = isDown & (g_start < b) & (g_end > a);
    if any(overlap)
        var(i) = 3;
    end
end
end

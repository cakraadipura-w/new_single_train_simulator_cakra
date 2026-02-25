% decision_var_NO_base.m
% Original (static) decision variable tagging for each speed-limit section.
% Output: global var (values in {1,2,3}), dimension = sum(var)

global var vel_profile gradient Davis Mass gravity

if isempty(gravity), gravity = 9.81; end
if isempty(var)
    var = zeros(1, size(vel_profile,1)-1);
end

% Reference speed for resistance estimate
v_ref = 30 / 3.6;
min_steepdownhill = -(Davis(1) + Davis(2)*v_ref + Davis(3)*v_ref^2) / Mass / gravity * 1000;

for i = 2:length(vel_profile)
    for j = 2:length(gradient)
        if gradient(j-1,2) <= min_steepdownhill && vel_profile(i-1,1) >= gradient(j-1,1) && vel_profile(i,1) <= gradient(j,1)
            var(i-1) = 1; % speed limit section fully inside steep downhill
            break;
        elseif gradient(j-1,2) <= min_steepdownhill && ( ...
                (vel_profile(i-1,1) <= gradient(j-1,1) && vel_profile(i,1) >= gradient(j-1,1)) || ...
                (vel_profile(i-1,1) <= gradient(j,1)   && vel_profile(i,1) >= gradient(j,1)) )
            var(i-1) = 3; % section overlaps steep downhill boundary
            break;
        else
            var(i-1) = 2; % no steep downhill
        end
    end
end

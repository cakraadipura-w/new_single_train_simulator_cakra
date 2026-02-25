% decision_var_NO.m  (DISPATCHER)
% Switch between base vs improved decision-variable tagging.
% Controlled by global use_improved.
%
% Required files on path:
%   - decision_var_NO_base.m
%   - decision_var_NO_improved.m

global use_improved vel_profile gradient var
global vel_profile_raw gradient_raw

% Default: base
if isempty(use_improved)
    use_improved = false;
end

if ~logical(use_improved)
    % restore raw route definitions if we previously ran improved mode
    if ~isempty(vel_profile_raw)
        vel_profile = vel_profile_raw;
    end
    if ~isempty(gradient_raw)
        gradient = gradient_raw;
    end

    % IMPORTANT: ensure var length matches current vel_profile (may have been refined previously)
    var = []; % let base file size it correctly

    decision_var_NO_base;
else
    decision_var_NO_improved;

    % (optional safety) ensure var length matches refined vel_profile
    if ~isempty(vel_profile)
        nSec = size(vel_profile,1)-1;
        if numel(var) ~= nSec
            var = var(1:min(numel(var),nSec));
            if numel(var) < nSec
                var(end+1:nSec) = 2;
            end
        end
    end
end

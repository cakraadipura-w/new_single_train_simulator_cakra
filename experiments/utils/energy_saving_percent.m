function saving_percent = energy_saving_percent(E_proposed, E_baseline)
%ENERGY_SAVING_PERCENT  Percentage energy saving relative to flat-out baseline.
%   saving_percent = (E_baseline - E_proposed) / E_baseline * 100
% Positive = saving, negative = worse. Supports scalar or vector inputs.

    saving_percent = (E_baseline - E_proposed) ./ (E_baseline + eps) * 100;
end

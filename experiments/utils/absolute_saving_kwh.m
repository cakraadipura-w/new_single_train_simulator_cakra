function delta_kwh = absolute_saving_kwh(E_proposed, E_baseline)
%ABSOLUTE_SAVING_KWH  Absolute energy saving (kWh) vs flat-out baseline.
%   Positive = saving. Supports scalar or vector inputs.

    delta_kwh = E_baseline - E_proposed;
end

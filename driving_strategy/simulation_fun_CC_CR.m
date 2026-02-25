function [running_inter, Total_E, s_out, v_out, vlim_out] = simulation_fun_CC_CR(X)
    %SIMULATION_FUN_CC_CR Wrapper to switch base vs improved simulator.
    % Controlled by global use_improved.
    %
    % Required:
    %   - simulation_fun_CC_CR_base.m
    %   - simulation_fun_CC_CR_improved.m
    
    global use_improved
    
    if isempty(use_improved) || ~logical(use_improved)
        [running_inter, Total_E, s_out, v_out, vlim_out] = simulation_fun_CC_CR_base(X);
    else
        [running_inter, Total_E, s_out, v_out, vlim_out] = simulation_fun_CC_CR_improved(X);
    end
end

function [pop, GD, SP] = nsga2_main(vel_profile)
%NSGA2_MAIN DISPATCHER - Choose NSGA-II variant via global nsga2_variant.
%
% This dispatcher provides a unified interface for multiple NSGA-II variants:
%   - 'original'   : Standard NSGA-II (crowding distance)
%   - 'rl_sde'     : RL-based mate selection + Shift-Density Estimation
%   - 'improved_regular' : Improved NSGA-II (regular variant)
%   - 'bayes_rl'   : Bayesian Reinforcement Learning only
%
% All variants are normalized to return exactly 3 outputs: [pop, GD, SP]
%
% Configuration: Set global nsga2_variant or cfg.nsga2_variant in main.m

    global nsga2_variant

    % Default variant
    if isempty(nsga2_variant)
        nsga2_variant = 'rl_sde';
    end
    
    variant = lower(string(nsga2_variant));
    
    switch variant
        case 'original'
            % Standard NSGA-II with crowding distance
            % NOTE: nsga2_main_original returns [G × N × (D+4)] — squeeze to final generation
            [allpop, GD, SP] = nsga2_main_original(vel_profile);
            pop = squeeze(allpop(end, :, :));  % Extract last generation → [N × (D+4)]
            
        case 'rl_sde'
            % RL + SDE hybrid (extracts first 3 outputs, discards history)
            [pop_extended, GD, SP, ~] = nsga2_main_rl_sde(vel_profile);
            pop = pop_extended;
            
        case 'improved_regular'
            % Improved NSGA-II (regular variant)
            [pop, GD, SP] = nsga2_main_improved_regular(vel_profile);
            
        case 'bayes_rl'
            % Bayesian RL variant
            [pop, GD, SP] = nsga2_main_bayes_RL(vel_profile);
            
        otherwise
            error('Unknown NSGA-II variant: "%s"\nAvailable: original, rl_sde, improved_regular, bayes_rl', variant);
    end
end
     

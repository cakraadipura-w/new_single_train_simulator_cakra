function OUT = run_single(cfg)
    %RUN_SINGLE Run one optimizer using current globals.
    % Returns: OUT struct with optimizer result + metadata for better file naming
    
    global pop_size iterations vel_profile use_improved nsga2_variant
    pop_size   = cfg.pop_size;
    iterations = cfg.iterations;
    
    opt = lower(string(cfg.optimizer));
    disp("Running single optimizer: " + opt);
    
    t0 = tic;
    switch opt
        case "nsga2"
            [pop, GD, SP] = nsga2_main(vel_profile);
        case "mopso"
            pop = mopso_main(vel_profile); GD = []; SP = [];
        case "spea2"
            pop = spea2_main(vel_profile); GD = []; SP = [];
        case "moead"
            pop = moead_main(vel_profile); GD = []; SP = [];
        otherwise
            error('Unknown optimizer: %s', opt);
    end
    rt = toc(t0);
    
    % Determine variant and strategy for result naming
    if strcmp(opt, 'nsga2')
        variant = lower(string(nsga2_variant));
    else
        variant = 'default';
    end
    
    strategy = 'base';
    if ~isempty(use_improved) && logical(use_improved)
        strategy = 'improved';
    end
    
    OUT = struct();
    OUT.optimizer = opt;
    OUT.nsga2_variant = variant;  % Added for better identification
    OUT.strategy = strategy;        % Added for better identification
    OUT.runtime_s = rt;
    OUT.pop = pop;
    OUT.GD  = GD;
    OUT.SP  = SP;
    
    fprintf('%s (%s) [%s] done | runtime=%.2fs\n', opt, variant, strategy, rt);
end

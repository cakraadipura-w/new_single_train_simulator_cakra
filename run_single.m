function OUT = run_single(cfg)
    %RUN_SINGLE Run one optimizer using current globals.
    
    global pop_size iterations vel_profile
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
    
    OUT = struct();
    OUT.optimizer = opt;
    OUT.runtime_s = rt;
    OUT.pop = pop;
    OUT.GD  = GD;
    OUT.SP  = SP;
    
    fprintf('%s done | runtime=%.2fs\n', opt, rt);
end

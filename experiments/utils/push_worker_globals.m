function push_worker_globals(ps, it, tobj, nsga2_var)
%PUSH_WORKER_GLOBALS  Set optimisation globals on each parallel worker.
%
% Usage:
%   wait(parfevalOnAll(@push_worker_globals, 0, POP_SIZE, ITERATIONS, T_TARGET, 'rl_sde'));
%
% Called via parfevalOnAll so it runs inside each worker's MATLAB session.

    global pop_size iterations time_obj_max nsga2_variant
    pop_size   = ps;
    iterations = it;
    if nargin >= 3 && ~isempty(tobj)
        time_obj_max = tobj;
    end
    if nargin >= 4 && ~isempty(nsga2_var)
        nsga2_variant = nsga2_var;
    end
end

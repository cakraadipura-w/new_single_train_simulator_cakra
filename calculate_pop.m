function evaluate = calculate_pop(pop)
%CALCULATE_POP Unified evaluator for any MOEA solver.
% pop: N x dim  -> evaluate: N x 2 [T, E]

global parallel_use driving_strategy
global show_progress


% pilih simulator sekali (aman buat parfor)
simfun = pick_simfun(driving_strategy);

n   = size(pop,1);
dim = size(pop,2);

fx1 = zeros(n,1);
fx2 = zeros(n,1);

pool   = gcp('nocreate');
usePar = logical(parallel_use) && ~isempty(pool);

doLog = false;
if exist('show_progress','var') && ~isempty(show_progress)
    doLog = logical(show_progress);
end

if usePar
    if doLog
        dq = parallel.pool.DataQueue;
        t0 = tic;

        % ✅ pakai anonymous function (NO local function di tengah file)
        %afterEach(dq, @(msg) fprintf('eval %d/%d | elapsed %.1fs | worker %d\n',msg(1), n, toc(t0), msg(2)));

        parfor i = 1:n
            [fx1(i), fx2(i)] = simfun(pop(i,1:dim));
            % --- sanitize outputs (avoid stuck/complex/NaN solutions dominating) ---
            t = fx1(i); e = fx2(i);
            if ~isreal(t) || ~isreal(e); t = real(t); e = real(e); end
            if ~isfinite(t) || ~isfinite(e) || t<=0 || e<0 || t>1000
                t = 1e6; e = 1e6;
            end
            fx1(i)=t; fx2(i)=e;

            if mod(i,20)==0 || i==n
                task = getCurrentTask();
                if isempty(task), wid = 0; else, wid = task.ID; end
                send(dq, [i, wid]);
            end
        end
    else
        parfor i = 1:n
            [fx1(i), fx2(i)] = simfun(pop(i,1:dim));
            % --- sanitize outputs (avoid stuck/complex/NaN solutions dominating) ---
            t = fx1(i); e = fx2(i);
            if ~isreal(t) || ~isreal(e); t = real(t); e = real(e); end
            if ~isfinite(t) || ~isfinite(e) || t<=0 || e<0 || t>1000
                t = 1e6; e = 1e6;
            end
            fx1(i)=t; fx2(i)=e;
        end
    end
else
    t0 = tic;
    for i = 1:n
        [fx1(i), fx2(i)] = simfun(pop(i,1:dim));
        % --- sanitize outputs (avoid stuck/complex/NaN solutions dominating) ---
        t = fx1(i); e = fx2(i);
        if ~isreal(t) || ~isreal(e); t = real(t); e = real(e); end
        if ~isfinite(t) || ~isfinite(e) || t<=0 || e<0 || t>1000
            t = 1e6; e = 1e6;
        end
        fx1(i)=t; fx2(i)=e;
        if doLog && (mod(i,20)==0 || i==n)
            %fprintf('eval %d/%d | elapsed %.1fs (serial)\n', i, n, toc(t0));
        end
    end
end

evaluate = [fx1, fx2];
end

function simfun = pick_simfun(ds)
    ds = upper(string(ds));
    global use_improved
    switch ds
        case "CC"
            simfun = @simulation_fun_CC;
        case "CR"
            simfun = @simulation_fun_CR;
        case "CC_CR"
             simfun = @simulation_fun_CC_CR;
        otherwise
            error('unknown driving_strategy: %s', ds);
    end
end
function [pop, GD, SP] = de_moea_main(vel_profile)
%DE_MOEA_MAIN  Multi-Objective Differential Evolution (MODE) with NSGA-II selection.
%
% This is the "Improved DE" variant referenced as [N13] in the experiment spec (E3).
% Combines standard DE mutation/crossover operators with NSGA-II non-dominated
% sorting and crowding distance selection.
%
% DE parameters:
%   F_scale = 0.5  (differential weight)
%   CR      = 0.9  (crossover probability)
%   Variant : DE/rand/1/bin
%
% Returns:
%   pop  : N × (dim+4) [X | T | E | rank | crowd]  (same format as NSGA-II)
%   GD   : 1×G generalised distance (NaN — not tracked)
%   SP   : 1×G spacing (NaN — not tracked)
%
% Requires globals: pop_size, iterations, dimension, var

    global pop_size iterations dimension var
    global parallel_use

    if isempty(pop_size) || isempty(iterations) || isempty(dimension)
        error('de_moea_main: globals pop_size/iterations/dimension not set.');
    end

    N  = pop_size;
    G  = iterations;

    % --- DE hyper-parameters ---
    F_scale = 0.5;   % differential weight
    CR      = 0.9;   % crossover probability

    % --- logging / display ---
    log_each    = true;
    global show_progress
    plotDisplay = ~isempty(show_progress) && isscalar(show_progress) && logical(show_progress);
    plotEvery   = 20;

    % --- bounds (same as mopso_main) ---
    lb = 20 + zeros(1, dimension);
    ub = lb;
    for i = 1:dimension
        for j = 1:length(var)
            if j == 1 && i <= var(j)
                ub(i) = vel_profile(j, 2);
            elseif i <= sum(var(1:j)) && i > sum(var(1:j-1))
                ub(i) = vel_profile(j, 2);
                break;
            end
        end
    end

    % --- initialise population ---
    popX = lb + rand(N, dimension) .* (ub - lb);
    popF = calculate_pop(popX);

    % --- initial NSGA-II ranking ---
    [popRank, popCrowd] = rank_and_crowd(popF);

    GD = nan(1, G);
    SP = nan(1, G);

    tStart = tic;

    for gen = 1:G
        % =================================================================
        % DE/rand/1/bin  mutation + crossover
        % =================================================================
        childX = zeros(N, dimension);
        for i = 1:N
            % Choose 3 distinct random individuals (all different from i)
            candidates = setdiff(1:N, i);
            r = candidates(randperm(length(candidates), 3));
            r1 = r(1); r2 = r(2); r3 = r(3);

            % Mutation: V = X_r1 + F*(X_r2 - X_r3)
            V = popX(r1,:) + F_scale .* (popX(r2,:) - popX(r3,:));
            V = max(lb, min(ub, V));  % clamp to bounds

            % Crossover: bin (binomial)
            jrand = randi(dimension);  % at least one gene comes from V
            mask = rand(1, dimension) < CR;
            mask(jrand) = true;
            U = popX(i,:);
            U(mask) = V(mask);

            childX(i,:) = U;
        end

        % --- evaluate offspring ---
        childF = calculate_pop(childX);

        % --- NSGA-II merge + select ---
        combX = [popX; childX];
        combF = [popF;  childF];
        [combRank, combCrowd] = rank_and_crowd(combF);
        [popX, popF, popRank, popCrowd] = select_next_gen(combX, combF, combRank, combCrowd, N);

        % --- logging ---
        if log_each && (mod(gen,20)==0 || gen==1 || gen==G)
            elapsed = toc(tStart);
            front1  = popF(popRank == 1, :);
            nF1     = size(front1, 1);
            [bestE, iE] = min(front1(:,2)); T_at_E = front1(iE,1);
            [bestT, iT] = min(front1(:,1)); E_at_T = front1(iT,2);
            fprintf('[DE-MOEA] gen %4d/%d | F1=%d | bestE=%.4f (T=%.1f s) | bestT=%.1f (E=%.4f) | %.1f s\n', ...
                gen, G, nF1, bestE, T_at_E, bestT, E_at_T, elapsed);
        end

        % --- plot ---
        if plotDisplay && (mod(gen, plotEvery)==0 || gen==1)
            figure(201); clf;
            xP = real(popF(:,1)); yP = real(popF(:,2));
            mP = isfinite(xP) & isfinite(yP);
            scatter(xP(mP), yP(mP), 'filled');
            grid on;
            xlabel('Running time (s)'); ylabel('Energy (kWh)');
            title(sprintf('DE-MOEA gen %d', gen));
            drawnow;
        end
    end

    % --- assemble output (same format as NSGA-II) ---
    pop = [popX, popF, popRank, popCrowd];
end

% =========================================================================
%  LOCAL HELPERS (identical interface to nsga2_main_original helpers)
% =========================================================================

function [rank, crowd] = rank_and_crowd(F)
    n = size(F,1);
    rank  = zeros(n,1);
    crowd = zeros(n,1);

    remaining = 1:n;
    r = 1;
    while ~isempty(remaining)
        nd = nondom_indices(F(remaining,:));
        front_idx = remaining(nd);
        rank(front_idx) = r;
        crowd(front_idx) = crowding_dist(F(front_idx,:));
        remaining = remaining(~nd);
        r = r + 1;
    end
end

function nd = nondom_indices(F)
    n  = size(F,1);
    nd = true(n,1);
    for i = 1:n
        if ~nd(i), continue; end
        for j = 1:n
            if i==j || ~nd(j), continue; end
            if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                nd(i) = false; break;
            end
        end
    end
end

function cd = crowding_dist(F)
    n = size(F,1);
    if n <= 2, cd = Inf(n,1); return; end
    cd = zeros(n,1);
    for m = 1:size(F,2)
        [~, idx] = sort(F(:,m));
        cd(idx(1)) = Inf; cd(idx(end)) = Inf;
        frange = F(idx(end),m) - F(idx(1),m) + eps;
        for i = 2:n-1
            cd(idx(i)) = cd(idx(i)) + (F(idx(i+1),m) - F(idx(i-1),m)) / frange;
        end
    end
end

function [nX, nF, nR, nC] = select_next_gen(combX, combF, combRank, combCrowd, N)
    n = size(combF,1);
    order = zeros(n,1);
    ptr   = 1;
    r = 1;
    while ptr <= N
        idx_r = find(combRank == r);
        if ptr + length(idx_r) - 1 <= N
            order(ptr:ptr+length(idx_r)-1) = idx_r;
            ptr = ptr + length(idx_r);
        else
            need = N - ptr + 1;
            [~, srt] = sort(combCrowd(idx_r), 'descend');
            chosen = idx_r(srt(1:need));
            order(ptr:N) = chosen;
            ptr = N + 1;
        end
        r = r + 1;
    end
    sel = order(1:N);
    nX = combX(sel,:); nF = combF(sel,:);
    nR = combRank(sel); nC = combCrowd(sel);
end

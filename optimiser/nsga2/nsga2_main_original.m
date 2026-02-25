function [allpop, GD, SP] = nsga2_main_original(vel_profile)
%NSGA2_MAIN One-file NSGA-II (minimize [T, E]) calling external calculate_pop(popX).
% Keeps compatibility with run_single.m calling nsga2_main(vel_profile).
%
% Requires globals set in setup_project/setup_globals_parallel:
%   pop_size, iterations, dimension
%
% Evaluator:
%   calculate_pop(popX) -> N x 2 [T, E]
%
% Notes:
% - Bounds are set generically from vel_profile max speed limit (km/h) and a floor (20 km/h).
%   This avoids mismatch between gene->section mapping while remaining safe for optimizers.
% - Non-dominated sorting uses correct minimization dominance:
%     A dominates B if all(A<=B) && any(A<B)

    %#ok<*NASGU>
    global pop_size iterations dimension

    if isempty(pop_size) || isempty(iterations) || isempty(dimension)
        error('Globals pop_size/iterations/dimension not set. Check setup_project.m.');
    end

    N  = pop_size;
    G  = iterations;
    M  = 2; % objectives

    % --- user toggles
    log_each = true;
    plotDisplay = true;
    plotEvery = 10;

    % --- bounds (km/h)
    vmax = max(vel_profile(:,2));
    lb = 20 * ones(1, dimension);
    ub = vmax * ones(1, dimension);
    bounds = [lb(:), ub(:)];

    % --- init population
    popX = rand_init(N, dimension, bounds);
    popF = calculate_pop(popX);

    % --- rank+crowd
    [popRank, popCrowd] = rank_and_crowd(popF);

    % --- history storage: [X | T | E | rank | crowd]
    allpop = zeros(G, N, dimension+4);
    GD = zeros(1,G);
    SP = zeros(1,G);

    tStart = tic;

    for gen = 1:G
        % --- create offspring
        childX = make_children(popX, popRank, popCrowd, bounds);

        % --- evaluate offspring
        childF = calculate_pop(childX);

        % --- merge
        combX = [popX; childX];
        combF = [popF; childF];

        % --- sort merged
        [combRank, combCrowd] = rank_and_crowd(combF);

        % --- select next generation
        [popX, popF, popRank, popCrowd] = select_next_gen(combX, combF, combRank, combCrowd, N);

        % --- save
        allpop(gen,:,:) = [popX, popF, popRank, popCrowd];

        % --- SP metric (reuse your calculate_sp if it exists; else NaN)
        if exist('calculate_sp','file')==2
            try
                SP(gen) = calculate_sp([popX, popF, popRank, popCrowd]);
            catch
                SP(gen) = NaN;
            end
        else
            SP(gen) = NaN;
        end

        % --- logging
        if log_each && (mod(gen,10)==0 || gen==1 || gen==G)
            elapsed = toc(tStart);
            eta = elapsed/gen * (G-gen);

            front1 = popF(popRank==1,:);
            nF1 = size(front1,1);

            [bestE, idxE] = min(front1(:,2));  T_at_bestE = front1(idxE,1);
            [bestT, idxT] = min(front1(:,1));  E_at_bestT = front1(idxT,2);

            fprintf('gen %4d/%d | pop=%d | F1=%d | bestE=%.4f kWh (T=%.2f s) | bestT=%.2f s (E=%.4f kWh) | elapsed=%.1fs | ETA=%.1fs\n', ...
                gen, G, N, nF1, bestE, T_at_bestE, bestT, E_at_bestT, elapsed, eta);
        end

        % --- plot
        if plotDisplay && (mod(gen,plotEvery)==0 || gen==1)
            figure(101); clf;
            xP = real(popF(:,1)); yP = real(popF(:,2));
        mP = isfinite(xP) & isfinite(yP);
        scatter(xP(mP), yP(mP), 'filled');
            grid on
            xlabel('Running time (s)'); ylabel('Energy (kWh)');
            title(sprintf('NSGA-II original gen %d', gen));
            drawnow;
        end
    end
end

% ===================== helpers (subfunctions) =====================

function X = rand_init(N, D, bounds)
    lb = bounds(:,1)'; ub = bounds(:,2)';
    X = lb + rand(N,D).*(ub-lb);
end

function children = make_children(popX, rank, crowd, bounds)
    % Parameters (sensible defaults)
    N = size(popX,1);
    D = size(popX,2);
    eta_c = 15;      % SBX index
    eta_m = 20;      % poly mutation index
    p_c   = 0.9;     % crossover probability
    p_m   = 1/D;     % mutation probability per gene

    parentN = N; % produce N children
    children = zeros(parentN, D);

    for i=1:2:parentN
        p1 = tournament_select(rank, crowd);
        p2 = tournament_select(rank, crowd);

        x1 = popX(p1,:); x2 = popX(p2,:);

        if rand < p_c
            [c1, c2] = sbx_pair(x1, x2, bounds, eta_c);
        else
            c1 = x1; c2 = x2;
        end

        c1 = poly_mutate(c1, bounds, eta_m, p_m);
        c2 = poly_mutate(c2, bounds, eta_m, p_m);

        children(i,:) = c1;
        if i+1 <= parentN
            children(i+1,:) = c2;
        end
    end
end

function idx = tournament_select(rank, crowd)
    % binary tournament: prefer lower rank, then higher crowding
    n = numel(rank);
    a = randi(n); b = randi(n);
    if rank(a) < rank(b)
        idx = a;
    elseif rank(b) < rank(a)
        idx = b;
    else
        if crowd(a) > crowd(b)
            idx = a;
        else
            idx = b;
        end
    end
end

function [c1, c2] = sbx_pair(p1, p2, bounds, eta_c)
    D = numel(p1);
    lb = bounds(:,1)'; ub = bounds(:,2)';
    c1 = zeros(1,D); c2 = zeros(1,D);
    for j=1:D
        u = rand;
        if abs(p1(j)-p2(j)) < 1e-12
            c1(j) = p1(j); c2(j) = p2(j);
            continue;
        end
        if u <= 0.5
            beta = (2*u)^(1/(eta_c+1));
        else
            beta = (1/(2*(1-u)))^(1/(eta_c+1));
        end
        c1(j) = 0.5*((1+beta)*p1(j) + (1-beta)*p2(j));
        c2(j) = 0.5*((1-beta)*p1(j) + (1+beta)*p2(j));
        % clamp
        c1(j) = min(max(c1(j), lb(j)), ub(j));
        c2(j) = min(max(c2(j), lb(j)), ub(j));
    end
end

function x = poly_mutate(x, bounds, eta_m, p_m)
    D = numel(x);
    lb = bounds(:,1)'; ub = bounds(:,2)';
    for j=1:D
        if rand < p_m
            u = rand;
            if u < 0.5
                delta = (2*u)^(1/(eta_m+1)) - 1;
            else
                delta = 1 - (2*(1-u))^(1/(eta_m+1));
            end
            x(j) = x(j) + delta*(ub(j)-lb(j));
            x(j) = min(max(x(j), lb(j)), ub(j));
        end
    end
end

function [rank, crowd] = rank_and_crowd(F)
    % F: N x 2, minimize
    fronts = fast_nondominated_sort(F);
    N = size(F,1);
    rank = zeros(N,1);
    crowd = zeros(N,1);
    for fi = 1:numel(fronts)
        idx = fronts{fi};
        rank(idx) = fi;
        crowd(idx) = crowding_distance(F, idx);
    end
end

function fronts = fast_nondominated_sort(F)
    N = size(F,1);
    S = cell(N,1);
    n = zeros(N,1);
    fronts = {};
    F1 = [];

    for p=1:N
        Sp = [];
        np = 0;
        for q=1:N
            if p==q, continue; end
            if dominates(F(p,:), F(q,:))
                Sp(end+1) = q; %#ok<AGROW>
            elseif dominates(F(q,:), F(p,:))
                np = np + 1;
            end
        end
        S{p} = Sp;
        n(p) = np;
        if n(p)==0
            F1(end+1) = p; %#ok<AGROW>
        end
    end

    fronts{1} = F1;
    i = 1;
    while ~isempty(fronts{i})
        Q = [];
        for p = fronts{i}
            for q = S{p}
                n(q) = n(q) - 1;
                if n(q)==0
                    Q(end+1) = q; %#ok<AGROW>
                end
            end
        end
        i = i + 1;
        fronts{i} = Q;
    end
    % remove last empty
    if isempty(fronts{end})
        fronts(end) = [];
    end
end

function d = crowding_distance(F, idx)
    % Standard NSGA-II crowding, for minimization
    m = size(F,2);
    d = zeros(numel(idx),1);
    if numel(idx) <= 2
        d(:) = inf;
        return;
    end
    Fi = F(idx,:);
    for k=1:m
        [vals, order] = sort(Fi(:,k));
        d(order(1)) = inf;
        d(order(end)) = inf;
        vmin = vals(1); vmax = vals(end);
        if vmax - vmin < 1e-12
            continue;
        end
        for j=2:numel(idx)-1
            d(order(j)) = d(order(j)) + (vals(j+1)-vals(j-1))/(vmax-vmin);
        end
    end
end

function flag = dominates(a, b)
    flag = all(a <= b) && any(a < b);
end

function [X, F, rank, crowd] = select_next_gen(combX, combF, combRank, combCrowd, N)
    % Sort by rank asc, crowd desc
    idx = (1:size(combX,1))';
    T = table(idx, combRank, combCrowd);
    T = sortrows(T, {'combRank','combCrowd'}, {'ascend','descend'});
    sel = T.idx(1:N);

    X = combX(sel,:);
    F = combF(sel,:);
    rank = combRank(sel,:);
    crowd = combCrowd(sel,:);
end

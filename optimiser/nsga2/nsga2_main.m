function [pop, GD, SP, history] = nsga2_main(vel_profile)
% NSGA2_MAIN_HYBRID (RL + SDE Improved NSGA-II)
% Kombinasi 2 Paper:
% 1. Bayesian Reinforcement Learning (Wen et al., 2025) untuk Seleksi Pasangan (Kasta).
% 2. Shift-Based Density Estimation / SDE (Zhang et al., 2025) pengganti Crowding Distance.
% Fokus 2 Objektif: [Time, Energy]

    global pop_size iterations dimension

    if isempty(pop_size) || isempty(iterations) || isempty(dimension)
        error('Globals pop_size/iterations/dimension not set.');
    end

    N  = pop_size;
    G  = iterations;
    
    % --- User toggles
    log_each = true;
    plotDisplay = true;
    plotEvery = 10;

    % --- Bounds per decision var (km/h), based on local speed limits
    global var
    bounds = build_bounds_cccr(vel_profile, var, 5); % 5 km/h floor

    % =====================================================================
    % 1. INISIALISASI PROBABILITAS RL (Bayesian Prior) dari Paper 1
    % =====================================================================
    RL_Probs = (1/3) * ones(3, 3);

    % --- Init populasi
    popX = rand_init(N, dimension, bounds);
    popF = calculate_pop(popX);

    % --- Rank & SDE (Ganti Crowding Distance dengan SDE dari Paper 2)
    [popRank, popSDE] = rank_and_sde(popF);
    % --- Optional history storage (disable by default to save RAM)
    store_history = false;
    history = [];
    if store_history
        history = zeros(G, N, dimension+4);
    end
    GD = zeros(1,G);
    SP = zeros(1,G);

    tStart = tic;

    for gen = 1:G
        % =================================================================
        % 2. STRATIFIKASI POPULASI 
        % Urutkan berdasarkan Rank (Asc) lalu SDE (Desc -> nilai tinggi lebih baik)
        % =================================================================
        idx_all = (1:N)';
        sort_table = sortrows([idx_all, popRank, popSDE], [2, -3]);
        sorted_indices = sort_table(:,1);
        
        nE = floor(N/3);
        nI = floor(N/3);
        idx_E = sorted_indices(1 : nE);
        idx_I = sorted_indices(nE+1 : nE+nI);
        idx_L = sorted_indices(nE+nI+1 : end);
        
        tier_map = zeros(N, 1);
        tier_map(idx_E) = 1;
        tier_map(idx_I) = 2;
        tier_map(idx_L) = 3;

        % =================================================================
        % 3. MEMBUAT OFFSPRING DENGAN SMART MATE SELECTION (RL)
        % =================================================================
        [childX, parent_records] = make_children_rl(popX, popRank, popSDE, bounds, tier_map, RL_Probs, idx_E, idx_I, idx_L);

        % --- Evaluasi offspring
        childF = calculate_pop(childX);

        % =================================================================
        % 4. EVALUASI HASIL CROSSOVER & BAYESIAN UPDATE (Paper 1)
        % =================================================================
        success_matrix = zeros(3,3); 
        
        for k = 1:size(parent_records, 1)
            t1 = parent_records(k, 1); 
            t2 = parent_records(k, 2); 
            p1_idx = parent_records(k, 3);
            p2_idx = parent_records(k, 4);
            
            p1_F = popF(p1_idx, :);
            p2_F = popF(p2_idx, :);
            c1_F = childF(2*k-1, :);
            
            if 2*k <= N
                c2_F = childF(2*k, :);
            else
                c2_F = [inf, inf]; 
            end
            
            % Expected Pairing: Anak mendominasi KEDUA orang tuanya
            is_c1_expected = dominates(c1_F, p1_F) && dominates(c1_F, p2_F);
            is_c2_expected = dominates(c2_F, p1_F) && dominates(c2_F, p2_F);
            
            if is_c1_expected || is_c2_expected
                success_matrix(t1, t2) = success_matrix(t1, t2) + 1;
            end
        end
        
        % Bayesian Update
        epsilon = 0.01; 
        for i = 1:3
            denom = 0;
            for j = 1:3
                expected_count = success_matrix(i,j) + epsilon;
                denom = denom + (RL_Probs(i,j) * expected_count);
            end
            for j = 1:3
                expected_count = success_matrix(i,j) + epsilon;
                RL_Probs(i,j) = (RL_Probs(i,j) * expected_count) / denom;
            end
            RL_Probs(i,:) = max(RL_Probs(i,:), 0.05);
            RL_Probs(i,:) = RL_Probs(i,:) / sum(RL_Probs(i,:)); 
        end

        % =================================================================
        % 5. ENVIRONMENTAL SELECTION DENGAN SDE
        % =================================================================
        combX = [popX; childX];
        combF = [popF; childF];

        [combRank, combSDE] = rank_and_sde(combF);
        [popX, popF, popRank, popSDE] = select_next_gen(combX, combF, combRank, combSDE, N);
        % --- Save to history (optional)
        if store_history
            history(gen,:,:) = [popX, popF, popRank, popSDE];
        end

        if exist('calculate_sp','file')==2
            try SP(gen) = calculate_sp([popX, popF, popRank, popSDE]); catch, SP(gen) = NaN; end
        else
            SP(gen) = NaN;
        end

        % --- Logging
        if log_each && (mod(gen,10)==0 || gen==1 || gen==G)
            elapsed = toc(tStart);
            eta = elapsed/gen * (G-gen);
            front1 = popF(popRank==1,:);
            nF1 = size(front1,1);
            [bestE, idxE] = min(front1(:,2));  T_at_bestE = front1(idxE,1);
            [bestT, idxT] = min(front1(:,1));  E_at_bestT = front1(idxT,2);

            fprintf('gen %4d/%d | pop=%d | F1=%d | bestE=%.4f (T=%.2f) | bestT=%.2f (E=%.4f) | ETA=%.1fs\n', ...
                gen, G, N, nF1, bestE, T_at_bestE, bestT, E_at_bestT, eta);
        end

        % --- Plot
        if plotDisplay && (mod(gen,plotEvery)==0 || gen==1)
            figure(101); clf;
            xP = real(popF(:,1)); yP = real(popF(:,2));
            mP = isfinite(xP) & isfinite(yP);
            scatter(xP(mP), yP(mP), 'filled');
            grid on
            xlabel('Running time (s)'); ylabel('Energy (kWh)');
            title(sprintf('Hybrid RL-SDE NSGA-II gen %d', gen));
            drawnow;
        end
    end

    % --- Return final population in standard format (N x (D+4))
    pop = [popX, popF, popRank, popSDE];
end

% ===================== SUBFUNCTIONS =====================

function [children, parent_records] = make_children_rl(popX, rank, sde_val, bounds, tier_map, RL_Probs, idx_E, idx_I, idx_L)
    N = size(popX,1);
    D = size(popX,2);
    eta_c = 15;      
    eta_m = 20;      
    p_c   = 0.9;     
    p_m   = 1/D;     

    children = zeros(N, D);
    parent_records = zeros(ceil(N/2), 4); 
    
    pair_count = 1;
    for i=1:2:N
        % Turnamen sekarang menggunakan SDE, bukan Crowding
        p1 = tournament_select(rank, sde_val);
        t1 = tier_map(p1); 
        
        r = rand();
        probs = RL_Probs(t1, :);
        cum_probs = cumsum(probs);
        if r <= cum_probs(1)
            t2 = 1; pool = idx_E;
        elseif r <= cum_probs(2)
            t2 = 2; pool = idx_I;
        else
            t2 = 3; pool = idx_L;
        end
        
        if isempty(pool)
            p2 = tournament_select(rank, sde_val);
            t2 = tier_map(p2);
        else
            p2 = pool(randi(length(pool))); 
        end

        x1 = popX(p1,:); x2 = popX(p2,:);

        if rand < p_c
            [c1, c2] = sbx_pair(x1, x2, bounds, eta_c);
        else
            c1 = x1; c2 = x2;
        end

        c1 = poly_mutate(c1, bounds, eta_m, p_m);
        c2 = poly_mutate(c2, bounds, eta_m, p_m);


        % Repair to enforce CCCR ordering constraints per section
        c1 = repair_cccr(c1, bounds);
        c2 = repair_cccr(c2, bounds);
        children(i,:) = c1;
        if i+1 <= N
            children(i+1,:) = c2;
        end
        
        parent_records(pair_count, :) = [t1, t2, p1, p2];
        pair_count = pair_count + 1;
    end
end

function X = rand_init(N, D, bounds)
    lb = bounds(:,1)'; ub = bounds(:,2)';
    X = zeros(N,D);
    for i = 1:N
        r = rand;
        if r < 0.55
            u = rand(1,D).^2;        % bias to low speeds (energy)
        elseif r < 0.85
            u = rand(1,D);           % uniform
        else
            u = 1 - rand(1,D).^2;    % bias to high speeds (time feasibility)
        end
        X(i,:) = lb + u.*(ub-lb);
        X(i,:) = repair_cccr(X(i,:), bounds);
    end
end

function idx = tournament_select(rank, metric_val)
    % metric_val = SDE (Lebih tinggi = lebih baik penyebarannya)
    n = numel(rank);
    a = randi(n); b = randi(n);
    if rank(a) < rank(b)
        idx = a;
    elseif rank(b) < rank(a)
        idx = b;
    else
        if metric_val(a) > metric_val(b)
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

function [rank, sde_val] = rank_and_sde(F)
    fronts = fast_nondominated_sort(F);
    N = size(F,1);
    rank = zeros(N,1);
    sde_val = zeros(N,1);
    for fi = 1:numel(fronts)
        idx = fronts{fi};
        rank(idx) = fi;
        sde_val(idx) = shift_density_estimation(F, idx);
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
                Sp(end+1) = q; 
            elseif dominates(F(q,:), F(p,:))
                np = np + 1;
            end
        end
        S{p} = Sp;
        n(p) = np;
        if n(p)==0
            F1(end+1) = p; 
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
                    Q(end+1) = q; 
                end
            end
        end
        i = i + 1;
        fronts{i} = Q;
    end
    if isempty(fronts{end})
        fronts(end) = [];
    end
end

function flag = dominates(a, b)
    flag = all(a <= b) && any(a < b);
end

% =====================================================================
% FUNGSI BARU: Shift-Based Density Estimation (Pengganti Crowding)
% Referensi: Paper 2 (Zhang et al., 2025)
% =====================================================================
function sde_val = shift_density_estimation(F, idx)
    k = length(idx);
    sde_val = zeros(k,1);
    if k <= 2
        sde_val(:) = inf; % Biarkan ujung-ujung pareto tetap bertahan
        return;
    end
    
    Fi = F(idx, :);
    
    % Normalisasi nilai fungsi objektif f'
    fmin = min(Fi, [], 1);
    fmax = max(Fi, [], 1);
    fmax(fmax == fmin) = fmax(fmax == fmin) + 1e-12; % Cegah Div-by-Zero
    F_norm = (Fi - fmin) ./ (fmax - fmin);
    
    % Hitung SDE untuk setiap individu di dalam Pareto Front ini
    for i = 1:k
        min_dist = inf;
        for j = 1:k
            if i == j, continue; end
            % Rumus shift: sde(f'_n(p_i), f'_n(p_j)) = max(0, f'_n(p_j) - f'_n(p_i))
            shift_diff = max(0, F_norm(j, :) - F_norm(i, :));
            dist = sqrt(sum(shift_diff.^2));
            
            if dist < min_dist
                min_dist = dist;
            end
        end
        sde_val(i) = min_dist;
    end
end

function [X, F, rank, sde_val] = select_next_gen(combX, combF, combRank, combSDE, N)
    idx = (1:size(combX,1))';
    T = table(idx, combRank, combSDE);
    % Sort by Rank (Ascending), lalu SDE (Descending)
    T = sortrows(T, {'combRank','combSDE'}, {'ascend','descend'});
    sel = T.idx(1:N);

    X = combX(sel,:);
    F = combF(sel,:);
    rank = combRank(sel,:);
    sde_val = combSDE(sel,:);
end


% =====================================================================
% CCCR-specific helpers: bounds + repair (keeps optimizer from wasting evals)
% =====================================================================
function bounds = build_bounds_cccr(vel_profile, var, lb_floor)
    if nargin < 3 || isempty(lb_floor), lb_floor = 5; end % km/h
    if isempty(var)
        vmax = max(vel_profile(:,2));
        lb = lb_floor * ones(1, size(vel_profile,1)-1);
        ub = vmax * ones(1, numel(lb));
        bounds = [lb(:), ub(:)];
        return;
    end

    nSec = numel(var);
    vlim = vel_profile(1:nSec, 2); % km/h per section (start row)
    D = sum(var);
    lb = zeros(D,1);
    ub = zeros(D,1);

    k = 1;
    for i = 1:nSec
        vmax_i = max(lb_floor, vlim(i));
        vi = var(i);
        if vi < 1, vi = 1; end
        for j = 1:vi
            lb(k) = lb_floor;
            ub(k) = vmax_i;
            k = k + 1;
        end
    end
    bounds = [lb, ub];
end

function x = repair_cccr(x, bounds)
    % Enforce ordering constraints:
    %  var=2: cruising >= coasting
    %  var=3: HIGH >= MID >= LOW1 (where MID is 3rd variable / cst_low2)
    global var
    lb = bounds(:,1)'; ub = bounds(:,2)';
    x = min(max(x, lb), ub);

    if isempty(var)
        return;
    end

    k = 1;
    for i = 1:numel(var)
        switch var(i)
            case 1
                k = k + 1;

            case 2
                cr = x(k); co = x(k+1);
                if co > cr
                    co = cr;
                end
                x(k) = cr; x(k+1) = co;
                k = k + 2;

            case 3
                hi  = x(k);
                lo1 = x(k+1);
                mid = x(k+2);

                if hi < lo1
                    tmp = hi; hi = lo1; lo1 = tmp;
                end
                mid = min(max(mid, lo1), hi);

                x(k)   = hi;
                x(k+1) = lo1;
                x(k+2) = mid;
                k = k + 3;

            otherwise
                % fallback: treat as 1 variable
                k = k + 1;
        end
    end

    % final clamp (repair might move mid)
    x = min(max(x, lb), ub);
end

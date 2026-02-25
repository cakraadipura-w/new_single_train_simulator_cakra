function [allpop, GD, SP] = nsga2_main_bayes_RL(vel_profile)
% NSGA2_MAIN_DRL (RL-Improved NSGA-II)
% Meminimalkan [T, E] dengan mekanisme Bayesian Reinforcement Learning 
% untuk strategi seleksi pasangan adaptif (Wen et al., 2025).
% 
% Pastikan globals: pop_size, iterations, dimension sudah diset.

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

    % --- Bounds (km/h) - Menyesuaikan dengan batas kecepatan
    vmax = max(vel_profile(:,2));
    lb = 20 * ones(1, dimension);
    ub = vmax * ones(1, dimension);
    bounds = [lb(:), ub(:)];

    % =====================================================================
    % 1. INISIALISASI PROBABILITAS RL (Bayesian Prior)
    % Matriks 3x3 untuk probabilitas seleksi antar kasta:
    % Baris = Kasta Parent 1 (1:Elite, 2:Intermediate, 3:Low)
    % Kolom = Kasta Parent 2 (1:Elite, 2:Intermediate, 3:Low)
    % Awalnya semua memiliki probabilitas 1/3 (seragam).
    % =====================================================================
    RL_Probs = (1/3) * ones(3, 3);

    % --- Init populasi
    popX = rand_init(N, dimension, bounds);
    popF = calculate_pop(popX);

    % --- Rank & Crowd
    [popRank, popCrowd] = rank_and_crowd(popF);

    % --- History storage
    allpop = zeros(G, N, dimension+4);
    GD = zeros(1,G);
    SP = zeros(1,G);

    tStart = tic;

    for gen = 1:G
        % =================================================================
        % 2. STRATIFIKASI POPULASI (Hierarchical Division)
        % Urutkan berdasarkan Rank (Asc) lalu Crowding (Desc)
        % =================================================================
        idx_all = (1:N)';
        sort_table = sortrows([idx_all, popRank, popCrowd], [2, -3]);
        sorted_indices = sort_table(:,1);
        
        % Bagi menjadi 3 kasta (Elite, Intermediate, Low-level)
        nE = floor(N/3);
        nI = floor(N/3);
        idx_E = sorted_indices(1 : nE);
        idx_I = sorted_indices(nE+1 : nE+nI);
        idx_L = sorted_indices(nE+nI+1 : end);
        
        % Mapping cepat untuk mengetahui kasta suatu individu
        tier_map = zeros(N, 1);
        tier_map(idx_E) = 1;
        tier_map(idx_I) = 2;
        tier_map(idx_L) = 3;

        % =================================================================
        % 3. MEMBUAT OFFSPRING DENGAN SMART MATE SELECTION (RL)
        % =================================================================
        [childX, parent_records] = make_children_rl(popX, popRank, popCrowd, bounds, tier_map, RL_Probs, idx_E, idx_I, idx_L);

        % --- Evaluasi offspring
        childF = calculate_pop(childX);

        % =================================================================
        % 4. EVALUASI HASIL CROSSOVER & BAYESIAN UPDATE
        % Mengecek "Expected Pairing" (Apakah anak mendominasi kedua orang tua?)
        % =================================================================
        success_matrix = zeros(3,3); % Mencatat jumlah sukses tiap kombinasi kasta
        
        for k = 1:size(parent_records, 1)
            t1 = parent_records(k, 1); % Kasta Parent 1
            t2 = parent_records(k, 2); % Kasta Parent 2
            p1_idx = parent_records(k, 3);
            p2_idx = parent_records(k, 4);
            
            p1_F = popF(p1_idx, :);
            p2_F = popF(p2_idx, :);
            c1_F = childF(2*k-1, :);
            
            % Anak 2 ada jika N genap
            if 2*k <= N
                c2_F = childF(2*k, :);
            else
                c2_F = [inf, inf]; 
            end
            
            % Syarat "Expected Pairing": Anak mendominasi KEDUA orang tuanya
            is_c1_expected = dominates(c1_F, p1_F) && dominates(c1_F, p2_F);
            is_c2_expected = dominates(c2_F, p1_F) && dominates(c2_F, p2_F);
            
            if is_c1_expected || is_c2_expected
                success_matrix(t1, t2) = success_matrix(t1, t2) + 1;
            end
        end
        
        % --- Bayesian Formula Update (Persamaan 24 dari paper diadaptasi)
        epsilon = 0.01; % Mencegah probabilitas runtuh ke 0 jika tidak ada sukses
        for i = 1:3
            denom = 0;
            for j = 1:3
                expected_count = success_matrix(i,j) + epsilon;
                denom = denom + (RL_Probs(i,j) * expected_count);
            end
            
            % Update probabilitas posterior
            for j = 1:3
                expected_count = success_matrix(i,j) + epsilon;
                RL_Probs(i,j) = (RL_Probs(i,j) * expected_count) / denom;
            end
            
            % Batasan minimal eksplorasi agar tidak stagnan di satu strategi
            RL_Probs(i,:) = max(RL_Probs(i,:), 0.05);
            RL_Probs(i,:) = RL_Probs(i,:) / sum(RL_Probs(i,:)); % Normalisasi
        end

        % =================================================================
        % 5. STANDAR NSGA-II: MERGE & SELECT
        % =================================================================
        combX = [popX; childX];
        combF = [popF; childF];

        [combRank, combCrowd] = rank_and_crowd(combF);
        [popX, popF, popRank, popCrowd] = select_next_gen(combX, combF, combRank, combCrowd, N);

        % --- Save ke memori
        allpop(gen,:,:) = [popX, popF, popRank, popCrowd];

        % --- Perhitungan metrik SP (Jika ada fungsinya)
        if exist('calculate_sp','file')==2
            try
                SP(gen) = calculate_sp([popX, popF, popRank, popCrowd]);
            catch
                SP(gen) = NaN;
            end
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
            title(sprintf('RL-NSGA-II gen %d', gen));
            drawnow;
        end
    end
end

% ===================== SUBFUNCTIONS (RL & NSGA-II) =====================

function [children, parent_records] = make_children_rl(popX, rank, crowd, bounds, tier_map, RL_Probs, idx_E, idx_I, idx_L)
    N = size(popX,1);
    D = size(popX,2);
    eta_c = 15;      
    eta_m = 20;      
    p_c   = 0.9;     
    p_m   = 1/D;     

    children = zeros(N, D);
    parent_records = zeros(ceil(N/2), 4); % Format: [Kasta_P1, Kasta_P2, Idx_P1, Idx_P2]
    
    pair_count = 1;
    for i=1:2:N
        % 1. Pilih Parent 1 (Turnamen mempertahankan elitism)
        p1 = tournament_select(rank, crowd);
        t1 = tier_map(p1); % Tipe kasta parent 1
        
        % 2. Pilih Tipe Parent 2 dengan Roulette Wheel berdasarkan RL_Probs
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
        
        % Fallback jika pool kosong (jarang terjadi tapi butuh pengaman)
        if isempty(pool)
            p2 = tournament_select(rank, crowd);
            t2 = tier_map(p2);
        else
            p2 = pool(randi(length(pool))); % Ambil acak dari kasta target
        end

        x1 = popX(p1,:); 
        x2 = popX(p2,:);

        % Proses Crossover & Mutasi (Continuous SBX & Poly Mutation)
        if rand < p_c
            [c1, c2] = sbx_pair(x1, x2, bounds, eta_c);
        else
            c1 = x1; c2 = x2;
        end

        c1 = poly_mutate(c1, bounds, eta_m, p_m);
        c2 = poly_mutate(c2, bounds, eta_m, p_m);

        children(i,:) = c1;
        if i+1 <= N
            children(i+1,:) = c2;
        end
        
        % Simpan rekaman pasangan untuk perhitungan Evaluasi Crossover di iterasi berikutnya
        parent_records(pair_count, :) = [t1, t2, p1, p2];
        pair_count = pair_count + 1;
    end
end

% === FUNGSI ASLI NSGA-II (TIDAK BERUBAH) ===

function X = rand_init(N, D, bounds)
    lb = bounds(:,1)'; ub = bounds(:,2)';
    X = lb + rand(N,D).*(ub-lb);
end

function idx = tournament_select(rank, crowd)
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

function d = crowding_distance(F, idx)
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
    idx = (1:size(combX,1))';
    T = table(idx, combRank, combCrowd);
    T = sortrows(T, {'combRank','combCrowd'}, {'ascend','descend'});
    sel = T.idx(1:N);

    X = combX(sel,:);
    F = combF(sel,:);
    rank = combRank(sel,:);
    crowd = combCrowd(sel,:);
end
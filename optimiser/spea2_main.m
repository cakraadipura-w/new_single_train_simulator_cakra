function [archive, hist] = spea2_main(vel_profile)
    % SPEA2 for 2 objectives: minimize [Time, Energy]
    % Output archive format: [X, T, E, rank, crowd] (rank/crowd from sort_pop for plotting)
    
    global pop_size iterations dimension var
    
    N = pop_size;
    G = iterations;
    archive_size = pop_size;
    target = 2;
    log = false;
    global show_progress
    plotDisplay = ~isempty(show_progress) && isscalar(show_progress) && logical(show_progress);
    
    % ----- bounds (same logic as nsga2_main) -----
    lb = 20 + zeros(1,dimension);
    ub = lb;
    for i=1:dimension
        for j=1:length(var)
            if j==1 && i<=var(j)
                ub(i)=vel_profile(j,2);
            elseif i<=sum(var(1:j)) && i>sum(var(1:j-1))
                ub(i)=vel_profile(j,2);
                break;
            end
        end
    end
    
    % ----- init population -----
    X = rand(N,dimension).*(ub-lb) + lb;
    F = calculate_pop(X);               % [T,E]
    A_X = zeros(0,dimension);           % archive decision
    A_F = zeros(0,2);                   % archive objectives
    
    hist.bestE = zeros(1,G);
    hist.bestT = zeros(1,G);
    hist.nF1   = zeros(1,G);
    
    tStart = tic;

    warning off all
    % ===== single figure (reuse) =====
    if plotDisplay
        figId = 103;
        figure(figId); clf;
        set(gcf,'Name','SPEA2 Live','NumberTitle','off');
        ax = gca; grid(ax,'on');
        xlabel(ax,'Running time (s)');
        ylabel(ax,'Energy (kWh)');
        hPlot = plot(ax, nan, nan, '*');   % handle plot
    end
    
    for gen = 1:G
        % 1) union U = P ∪ A
        U_X = [X; A_X];
        U_F = [F; A_F];
        U_n = size(U_X,1);
    
        % 2) SPEA2 fitness
        k = max(1, floor(sqrt(U_n)));   % common SPEA2 choice
        fit = spea2_fitness(U_F, k);    % lower is better
    
        % 3) Environmental selection -> new archive
        sel = find(fit < 1);
        if isempty(sel)
            % if none <1, keep best by fitness
            [~,idx] = sort(fit, 'ascend');
            sel = idx(1:min(archive_size, numel(idx)));
        end
    
        A_X = U_X(sel,:);
        A_F = U_F(sel,:);
        A_fit = fit(sel);
    
        if size(A_X,1) > archive_size
            [A_X, A_F] = spea2_truncate(A_X, A_F, archive_size);
        elseif size(A_X,1) < archive_size
            % fill archive with best remaining by fitness
            remain = setdiff((1:U_n)', sel);
            [~,ord] = sort(fit(remain), 'ascend');
            need = archive_size - size(A_X,1);
            pick = remain(ord(1:min(need, numel(ord))));
            A_X = [A_X; U_X(pick,:)];
            A_F = [A_F; U_F(pick,:)];
            A_fit = [A_fit; fit(pick)];
        end
    
        % 4) Mating selection from archive (binary tournament on fitness)
        parents_X = zeros(N, dimension);
        for i=1:N
            p1 = randi(size(A_X,1));
            p2 = randi(size(A_X,1));
            if A_fit(p1) < A_fit(p2)
                parents_X(i,:) = A_X(p1,:);
            else
                parents_X(i,:) = A_X(p2,:);
            end
        end
    
        % 5) Variation (SBX + polynomial mutation) -> offspring population
        childX = reproduce_sbx_pm(parents_X, lb, ub);
        X = childX;
        F = calculate_pop(X);
    
        % 6) Logging (use NSGA-II sort_pop only for rank/crowd & reporting)
        tmp = sort_pop([A_X, A_F], target, dimension);
        front1 = tmp(tmp(:,dimension+3)==1,:);
    
        if isempty(front1)
            front1 = tmp;
        end
    
        [bestE, idxE] = min(front1(:,dimension+2));
        hist.bestE(gen) = bestE;
        hist.bestT(gen) = min(front1(:,dimension+1));
        hist.nF1(gen)   = size(front1,1);
        

        if mod(gen,10)==0 || gen==1
            if log
                elapsed = toc(tStart);
                eta = elapsed/gen*(G-gen);
                T_at_bestE = front1(idxE, dimension+1);
                fprintf('SPEA2 gen %3d/%d | A=%d | F1=%d | bestE=%.4f kWh (T=%.2f s) | elapsed=%.1fs | ETA=%.1fs\n', ...
                    gen, G, size(A_X,1), size(front1,1), bestE, T_at_bestE, elapsed, eta);
            end
        
            if plotDisplay
                set(hPlot, 'XData', front1(:,dimension+1), 'YData', front1(:,dimension+2));
                title(ax, sprintf('SPEA2 Gen %d/%d', gen, G));
                drawnow limitrate
            end
        end


    end
    
    % return archive with rank/crowd (for plotting consistent with NSGA-II)
    archive = sort_pop([A_X, A_F], target, dimension);

end

% ====================== helpers inside same file ======================

function fit = spea2_fitness(F, k)
    % F: [n x 2] objectives (minimize)
    n = size(F,1);
    
    % dominance matrix
    dom = false(n,n);
    for i=1:n
        for j=1:n
            if i==j, continue; end
            dom(i,j) = all(F(i,:)<=F(j,:)) && any(F(i,:)<F(j,:));
        end
    end
    
    % strength S(i) = number dominated by i
    S = sum(dom,2);
    
    % raw fitness R(i) = sum of strengths of dominators of i
    R = zeros(n,1);
    for i=1:n
        dominators = find(dom(:,i));
        if ~isempty(dominators)
            R(i) = sum(S(dominators));
        end
    end
    
    % density D(i) using k-th nearest neighbor in objective space
    D = zeros(n,1);
    dist = squareform_pdist(F);
    dist = dist + diag(inf(n,1));
    sd = sort(dist,2,'ascend');
    sigma_k = sd(:, min(k, n-1));
    D = 1 ./ (sigma_k + 2);
    
    fit = R + D;
end

function [X2, F2] = spea2_truncate(X, F, maxN)
    % distance-based truncation in objective space (standard SPEA2 idea)
    n = size(F,1);
    dist = squareform_pdist(F);
    dist = dist + diag(inf(n,1));
    
    keep = true(n,1);
    while sum(keep) > maxN
        idx = find(keep);
        dsub = dist(idx,idx);
        % for each point, find nearest neighbor distance
        nn = min(dsub,[],2);
        [~,worstLocal] = min(nn);     % smallest distance -> densest -> remove
        keep(idx(worstLocal)) = false;
    end
    
    X2 = X(keep,:);
    F2 = F(keep,:);
end

function D = squareform_pdist(F)
    n = size(F,1);
    D = zeros(n,n);
    for i=1:n
        for j=i+1:n
            d = norm(F(i,:)-F(j,:));
            D(i,j)=d; D(j,i)=d;
        end
    end
end

function childX = reproduce_sbx_pm(parentsX, lb, ub)
    % SBX crossover + polynomial mutation (simple, robust)
    [N, dim] = size(parentsX);
    childX = zeros(N,dim);
    
    pc = 1.0;        % crossover probability
    pm = 1.0/dim;    % mutation per gene
    eta_c = 20;      % SBX index
    eta_m = 20;      % mutation index
    
    for i=1:2:N
        p1 = parentsX(i,:);
        if i+1<=N, p2 = parentsX(i+1,:); else, p2 = parentsX(randi(N),:); end
    
        c1 = p1; c2 = p2;
    
        if rand < pc
            for j=1:dim
                u = rand;
                if u <= 0.5
                    beta = (2*u)^(1/(eta_c+1));
                else
                    beta = (2*(1-u))^(-1/(eta_c+1));
                end
                c1(j) = 0.5*((1+beta)*p1(j) + (1-beta)*p2(j));
                c2(j) = 0.5*((1-beta)*p1(j) + (1+beta)*p2(j));
            end
        end
    
        % polynomial mutation
        c1 = poly_mutate(c1, lb, ub, pm, eta_m);
        c2 = poly_mutate(c2, lb, ub, pm, eta_m);
    
        childX(i,:) = min(max(c1,lb),ub);
        if i+1<=N
            childX(i+1,:) = min(max(c2,lb),ub);
        end
    end
end

function x = poly_mutate(x, lb, ub, pm, eta_m)
    dim = numel(x);
    for j=1:dim
        if rand < pm
            u = rand;
            if u < 0.5
                delta = (2*u)^(1/(eta_m+1)) - 1;
            else
                delta = 1 - (2*(1-u))^(1/(eta_m+1));
            end
            x(j) = x(j) + delta*(ub(j)-lb(j))*0.1; % scaled step
        end
    end
    x = min(max(x,lb),ub);
end

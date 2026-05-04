function archive = mopso_main(vel_profile)
    % MOPSO for 2 objectives: [Running time, Energy]
    % Output: archive = [X, T, E, rank, crowding] (same style as NSGA-II after sort_pop)
    
    global pop_size iterations dimension var

    log = true;
    global show_progress
    plotDisplay = ~isempty(show_progress) && isscalar(show_progress) && logical(show_progress);
    
    target = 2;
    
    % ---- bounds (sama seperti nsga2_main) ----
    lb = 20 + zeros(1,dimension);
    ub = lb;
    for i = 1:dimension
        for j = 1:length(var)
            if j==1 && i<=var(j)
                ub(i)=vel_profile(j,2);
            elseif i<=sum(var(1:j)) && i>sum(var(1:j-1))
                ub(i)=vel_profile(j,2);
                break;
            end
        end
    end
    
    % ---- MOPSO params ----
    swarm_size   = pop_size;     % biar fair sama NSGA-II
    max_iter     = iterations;
    
    w_max = 0.9;                 % inertia start
    w_min = 0.4;                 % inertia end
    c1 = 1.5;                    % cognitive
    c2 = 1.5;                    % social
    
    vmax = 0.2*(ub-lb);          % velocity clamp (aman)
    
    archive_size = pop_size;     % archive size (Pareto memory)
    
    % ---- init swarm ----
    X = rand(swarm_size,dimension).*(ub-lb) + lb;
    V = zeros(swarm_size,dimension);
    
    F = calculate_pop(X);        % [T, E] (ini sudah parallel kalau pool aktif)
    
    % pbest init
    pbestX = X;
    pbestF = F;
    
    % init archive from swarm
    pop0 = [X, F];
    pop0 = sort_pop(pop0, target, dimension);
    archive = truncate_by_rank_crowd(pop0, archive_size, dimension);
    
    % optional plot
    %figure();

    if plotDisplay
        figId = 102;                 % beda dari NSGA2 (101)
        figure(figId); clf;
        set(gcf,'Name','MOPSO Live','NumberTitle','off');
        ax = gca; grid(ax,'on');
        hPlot = plot(ax, nan, nan, '*');   % handle plot biar update tanpa bikin plot baru
        xlabel(ax,'Running time (s)');
        ylabel(ax,'Energy (kWh)');
    end
    
    tStart = tic;
    
    for it = 1:max_iter
        % inertia schedule
        w = w_max - (w_max-w_min)*(it-1)/(max_iter-1);
    
        % ---- choose leader for each particle from archive (bias crowding) ----
        leaderX = pick_leaders(archive, swarm_size, dimension);
    
        % ---- update velocity & position ----
        r1 = rand(swarm_size,dimension);
        r2 = rand(swarm_size,dimension);
    
        V = w*V + c1*r1.*(pbestX - X) + c2*r2.*(leaderX - X);
    
        % clamp velocity
        V = max(V, -vmax);
        V = min(V,  vmax);
    
        X = X + V;
    
        % clamp position
        X = max(X, lb);
        X = min(X, ub);
    
        % ---- evaluate ----
        F = calculate_pop(X);
    
        % ---- update pbest (dominance) ----
        for i = 1:swarm_size
            newDomOld = dominates(F(i,:), pbestF(i,:));
            oldDomNew = dominates(pbestF(i,:), F(i,:));
    
            if newDomOld || (~oldDomNew && rand < 0.5)
                pbestX(i,:) = X(i,:);
                pbestF(i,:) = F(i,:);
            end
        end
    
        % ---- update archive: archive U current swarm ----
        combo = [archive(:,1:dimension+2); X, F]; % drop rank/crowd dulu
        combo = sort_pop(combo, target, dimension);
        archive = truncate_by_rank_crowd(combo, archive_size, dimension);
    
        % ---- log + plot tiap 10 iter ----
        if mod(it,10)==0 || it==1

            if log == true
                front1 = archive(archive(:,dimension+3)==1,:);
                [bestE, idxE] = min(front1(:,dimension+2));
                T_at_bestE = front1(idxE, dimension+1);
        
                elapsed = toc(tStart);
                eta = elapsed/it*(max_iter-it);
        
                fprintf('MOPSO it %3d/%d | F1=%d | bestE=%.4f kWh (T=%.2f s) | elapsed=%.1fs | ETA=%.1fs\n', it, max_iter, size(front1,1), bestE, T_at_bestE, elapsed, eta);
            end
                if plotDisplay == true
                    cla;
                    plot(archive(:,dimension+1), archive(:,dimension+2), '*');
                    grid on
                    xlabel('Running time (s)');
                    ylabel('Energy (kWh)');
                    title(['MOPSO Gen ', num2str(it)]);
                    drawnow limitrate
                end
        end
    end  
end
    
    % ===================== helpers (inside same file) =====================
    
function tf = dominates(a, b)
    % Minimization dominance: a dominates b if a<=b in all and < in any
    tf = all(a <= b) && any(a < b);
end
    
function leaders = pick_leaders(archive, swarm_size, dimension)
    % pick from rank-1 preferred, weight by crowding distance
    
    front1 = archive(archive(:,dimension+3)==1,:);
    if isempty(front1)
        front1 = archive;
    end
    
    cd = front1(:,dimension+4);
    cd(~isfinite(cd)) = max(cd(isfinite(cd))+eps); % handle Inf
    cd(cd<0) = 0;
    
    if all(cd==0)
        prob = ones(size(cd))/numel(cd);
    else
        prob = cd/sum(cd);
    end
    
    leaders = zeros(swarm_size, dimension);
    for i=1:swarm_size
        idx = randsample(size(front1,1), 1, true, prob);
        leaders(i,:) = front1(idx,1:dimension);
    end
end
    
function pop_out = truncate_by_rank_crowd(pop_in, maxN, dimension)
    % sort by rank asc, crowding desc; truncate   
    sorted = sortrows(pop_in, [dimension+3, -(dimension+4)]);
    if size(sorted,1) > maxN
        pop_out = sorted(1:maxN,:);
    else
        pop_out = sorted;
    end
 end

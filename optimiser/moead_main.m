function [pop, front1, hist] = moead_main(vel_profile)
    % MOEA/D (Tchebycheff decomposition) for 2 objectives: [Time, Energy]
    % pop format: [X, T, E, rank, crowd] (rank/crowd from sort_pop)
    % front1 is subset of pop with rank=1
    
    global pop_size iterations dimension var
    
    N = pop_size;
    G = iterations;
    target = 2;

    log = false;
    plotDisplay = true;
    
    % ----- bounds (same as nsga2_main) -----
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
    
    % ----- weight vectors for 2 objectives -----
    w1 = linspace(0,1,N)';
    W = [w1, 1-w1];
    W(W==0) = 1e-6; % avoid exact zero for scalarization
    
    % neighborhood on weight vectors
    T = min(20, N);  % neighbor size
    D = squareform_pdist(W);
    [~, ord] = sort(D, 2, 'ascend');
    B = ord(:, 1:T);
    
    % ----- init population -----
    X = rand(N,dimension).*(ub-lb) + lb;
    F = calculate_pop(X);          % [T,E]
    z = min(F, [], 1);             % reference point
    
    % parameters
    delta = 0.9;                   % prob use neighborhood parents
    tStart = tic;
    
    hist.bestE = zeros(1,G);
    hist.bestT = zeros(1,G);
    hist.nF1   = zeros(1,G);

    warning off all
    if plotDisplay
        figId = 104;
        figure(figId); clf;
        set(gcf,'Name','MOEA/D Live','NumberTitle','off');
        ax = gca; grid(ax,'on');
        xlabel(ax,'Running time (s)');
        ylabel(ax,'Energy (kWh)');
        hPlot = plot(ax, nan, nan, '*');   % handle plot
    end
    
    for gen = 1:G
        for i = 1:N
            if rand < delta
                P = B(i,:);
            else
                P = 1:N;
            end
            p1 = P(randi(numel(P)));
            p2 = P(randi(numel(P)));
    
            % reproduction: SBX+PM to create one child
            y = sbx_one_child(X(p1,:), X(p2,:), lb, ub);
            f_y = calculate_pop(y);  % 1x2
    
            % update reference
            z = min(z, f_y);
    
            % update neighbors by Tchebycheff
            for jj = 1:numel(P)
                j = P(jj);
                if g_tcheby(f_y, W(j,:), z) <= g_tcheby(F(j,:), W(j,:), z)
                    X(j,:) = y;
                    F(j,:) = f_y;
                end
            end
        end
    
        % logging using rank/crowd from NSGA-II sort_pop (only for reporting)
        tmp = sort_pop([X, F], target, dimension);
        front = tmp(tmp(:,dimension+3)==1,:);
        if isempty(front), front = tmp; end
    
        [bestE, idxE] = min(front(:,dimension+2));
        hist.bestE(gen) = bestE;
        hist.bestT(gen) = min(front(:,dimension+1));
        hist.nF1(gen)   = size(front,1);
        
        
        if mod(gen,10)==0 || gen==1
            if log == true
                elapsed = toc(tStart);
                eta = elapsed/gen*(G-gen);
                T_at_bestE = front(idxE, dimension+1);
                fprintf('MOEA/D gen %3d/%d | F1=%d | bestE=%.4f kWh (T=%.2f s) | elapsed=%.1fs | ETA=%.1fs\n', gen, G, size(front,1), bestE, T_at_bestE, elapsed, eta);
            end
            if plotDisplay
                set(hPlot, 'XData', front(:,dimension+1), 'YData', front(:,dimension+2));
                title(ax, sprintf('MOEA/D Gen %d/%d', gen, G));
                drawnow limitrate
            end
        end
    end
    
    pop = sort_pop([X, F], target, dimension);
    front1 = pop(pop(:,dimension+3)==1, :);
    
    end
    
    % ====================== helpers inside same file ======================
    
    function val = g_tcheby(f, w, z)
    val = max(w .* abs(f - z));
    end
    
    function y = sbx_one_child(p1, p2, lb, ub)
    % one child from SBX + polynomial mutation
    dim = numel(p1);
    eta_c = 20;
    eta_m = 20;
    pm = 1.0/dim;
    
    y = p1;
    for j=1:dim
        u = rand;
        if u <= 0.5
            beta = (2*u)^(1/(eta_c+1));
        else
            beta = (2*(1-u))^(-1/(eta_c+1));
        end
        y(j) = 0.5*((1+beta)*p1(j) + (1-beta)*p2(j));
    end
    
    % mutation
    for j=1:dim
        if rand < pm
            u = rand;
            if u < 0.5
                delta = (2*u)^(1/(eta_m+1)) - 1;
            else
                delta = 1 - (2*(1-u))^(1/(eta_m+1));
            end
            y(j) = y(j) + delta*(ub(j)-lb(j))*0.1;
        end
    end
    
    y = min(max(y, lb), ub);
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

function pop_sorted = sort_pop(pop_in, target, dimension)
%SORT_POP Append Pareto rank and crowding, then sort by rank/crowding.
% Input:
%   pop_in    : [X | F] or [X | F | ...]
%   target    : number of objective columns in F
%   dimension : number of decision-variable columns in X
%
% Output:
%   pop_sorted: [X | F | rank | crowd], sorted by rank asc and crowd desc.

    if nargin < 3
        error('sort_pop requires pop_in, target, and dimension.');
    end

    required_cols = dimension + target;
    if size(pop_in, 2) < required_cols
        error('sort_pop expects at least %d columns, got %d.', required_cols, size(pop_in, 2));
    end

    if isempty(pop_in)
        pop_sorted = zeros(0, required_cols + 2);
        return;
    end

    base_pop = pop_in(:, 1:required_cols);
    F = base_pop(:, dimension+1:dimension+target);

    [rank, crowd] = rank_and_crowd(F);
    pop_with_meta = [base_pop, rank, crowd];
    pop_sorted = sortrows(pop_with_meta, [required_cols + 1, -(required_cols + 2)]);
end

function [rank, crowd] = rank_and_crowd(F)
    n = size(F, 1);
    rank = zeros(n, 1);
    crowd = zeros(n, 1);

    remaining = 1:n;
    front_rank = 1;
    while ~isempty(remaining)
        nd = nondom_indices(F(remaining, :));
        front_idx = remaining(nd);
        rank(front_idx) = front_rank;
        crowd(front_idx) = crowding_dist(F(front_idx, :));
        remaining = remaining(~nd);
        front_rank = front_rank + 1;
    end
end

function nd = nondom_indices(F)
    n = size(F, 1);
    nd = true(n, 1);

    for i = 1:n
        if ~nd(i)
            continue;
        end
        for j = 1:n
            if i == j || ~nd(j)
                continue;
            end
            if all(F(j, :) <= F(i, :)) && any(F(j, :) < F(i, :))
                nd(i) = false;
                break;
            end
        end
    end
end

function cd = crowding_dist(F)
    n = size(F, 1);
    if n <= 2
        cd = Inf(n, 1);
        return;
    end

    cd = zeros(n, 1);
    for m = 1:size(F, 2)
        [~, idx] = sort(F(:, m));
        cd(idx(1)) = Inf;
        cd(idx(end)) = Inf;

        frange = F(idx(end), m) - F(idx(1), m) + eps;
        for i = 2:n-1
            cd(idx(i)) = cd(idx(i)) + (F(idx(i+1), m) - F(idx(i-1), m)) / frange;
        end
    end
end
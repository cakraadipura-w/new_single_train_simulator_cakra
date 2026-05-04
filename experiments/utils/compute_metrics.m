function [hv, igd_val, spread, feas_rate] = compute_metrics(F_pop, F_ref, T_target)
%COMPUTE_METRICS  Multi-objective quality metrics for one run.
%
% Inputs:
%   F_pop    : N×2 [T(s), E(kWh)] — full final population (includes penalised)
%   F_ref    : M×2 [T, E] — reference Pareto front (union over all runs)
%              Pass [] to compute metrics relative to the run's own front.
%   T_target : scalar — feasibility threshold in seconds
%
% Outputs:
%   hv       : Hypervolume (higher = better)
%   igd_val  : Inverted Generational Distance (lower = better)
%   spread   : Delta spread, uniformity (lower = better)
%   feas_rate: Fraction of F_pop solutions with T <= T_target

    PENALTY = 1e5;

    % ---- feasibility rate ----
    valid    = F_pop(:,1) < PENALTY & F_pop(:,2) < PENALTY;
    feasible = valid & F_pop(:,1) <= T_target;
    feas_rate = sum(feasible) / size(F_pop, 1);

    % ---- Pareto front of valid solutions ----
    F_valid = F_pop(valid, :);
    if isempty(F_valid)
        hv = 0; igd_val = NaN; spread = NaN;
        return;
    end
    F_nd = nondom2d(F_valid);

    % ---- normalisation reference ----
    if ~isempty(F_ref) && size(F_ref,1) > 1
        all_pts = [F_nd; F_ref];
    else
        all_pts = F_nd;
        F_ref   = F_nd;
    end
    mins  = min(all_pts, [], 1);
    maxs  = max(all_pts, [], 1);
    rng_  = maxs - mins + eps;

    F_nd_n  = (F_nd  - mins) ./ rng_;
    F_ref_n = (F_ref - mins) ./ rng_;
    ref_pt  = ones(1,2) + 0.1;

    hv      = hv2d(F_nd_n, ref_pt);
    igd_val = igd_metric(F_ref_n, F_nd_n);
    spread  = spread_metric(F_nd_n);
end

% =========================================================================
function Fnd = nondom2d(F)
    n = size(F,1); keep = true(n,1);
    for i = 1:n
        if ~keep(i), continue; end
        for j = 1:n
            if i==j || ~keep(j), continue; end
            if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
                keep(i) = false; break;
            end
        end
    end
    Fnd = F(keep,:);
    [~,ord] = sort(Fnd(:,1),'ascend');
    Fnd = Fnd(ord,:);
end

function hv = hv2d(F, ref)
    if isempty(F), hv = 0; return; end
    F = nondom2d(F);
    F(:,2) = cummin(F(:,2));
    hv = 0; prevE = ref(2);
    for i = 1:size(F,1)
        if F(i,2) > prevE, continue; end
        hv = hv + (ref(1) - F(i,1)) * (prevE - F(i,2));
        prevE = F(i,2);
    end
end

function val = igd_metric(R, A)
    if isempty(R) || isempty(A), val = NaN; return; end
    d = zeros(size(R,1),1);
    for i = 1:size(R,1)
        dr = A - R(i,:);
        d(i) = min(sqrt(sum(dr.^2, 2)));
    end
    val = mean(d);
end

function sp = spread_metric(F)
    if size(F,1) < 2, sp = 0; return; end
    diffs  = diff(F, 1, 1);
    d      = sqrt(sum(diffs.^2, 2));
    d_mean = mean(d);
    if d_mean < eps, sp = 0; return; end
    sp = sum(abs(d - d_mean)) / (length(d) * d_mean);
end

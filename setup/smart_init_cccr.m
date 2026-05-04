function pop = smart_init_cccr(n_seeds, bounds)
% SMART_INIT_CCCR  Physics-informed initial population for improved CC_CR.
% Generates seeds where cruising >= coasting (var=2) and HIGH >= LOW1 (var=3).
%
% Strategies (25% each):
%   1. cruise-heavy  : high speed, mild coasting
%   2. coast-heavy   : lower speed, aggressive coasting (energy focus)
%   3. balanced      : moderate
%   4. perturbed     : noise around balanced (diversity)

    global var llm_advisor_enabled llm_advisor_profile

    dim = size(bounds, 1);
    lb  = bounds(:,1)';
    ub  = bounds(:,2)';
    rng_b = ub - lb;

    if isempty(var)
        pop = rand(n_seeds, dim) .* rng_b + lb;
        return;
    end

    if ~isempty(llm_advisor_enabled) && logical(llm_advisor_enabled)
        pop = advisor_seed_init(n_seeds, bounds, var, llm_advisor_profile);
        if ~isempty(pop)
            return;
        end
    end

    % [main_scale, coasting_ratio]
    % DV1 = lb + (ub-lb)*main_scale
    % DV2 = lb + (DV1-lb)*coasting_ratio
    strategies = [0.90, 0.75;
                  0.70, 0.55;
                  0.82, 0.68;
                  0.82, 0.68];
    noise_levels = [0.03, 0.03, 0.03, 0.08];

    n_per  = floor(n_seeds / 4);
    counts = [n_per, n_per, n_per, n_seeds - 3*n_per];

    pop = zeros(n_seeds, dim);
    row = 0;

    for s = 1:4
        sc = strategies(s, 1);
        rt = strategies(s, 2);
        ns = counts(s);
        nz = noise_levels(s);

        % build template
        tmpl = zeros(1, dim);
        k = 1;
        for i = 1:length(var)
            lo   = lb(k);
            hi_b = ub(k);
            dv1  = lo + (hi_b - lo) * sc;
            switch var(i)
                case 2
                    dv2 = lo + (dv1 - lo) * rt;
                    tmpl(k)   = dv1;
                    tmpl(k+1) = dv2;
                    k = k + 2;
                case 3
                    dv2 = lo + (dv1 - lo) * rt;
                    tmpl(k)   = dv1;
                    tmpl(k+1) = dv2;
                    k = k + 2;
                otherwise  % var=1 or unexpected
                    tmpl(k) = dv1;
                    k = k + 1;
            end
        end

        % replicate with noise
        for j = 1:ns
            x = tmpl + (rand(1, dim) - 0.5) .* rng_b * nz;
            x = min(max(x, lb), ub);
            % enforce DV1 >= DV2
            k = 1;
            for i = 1:length(var)
                switch var(i)
                    case {2, 3}
                        if x(k+1) > x(k), x(k+1) = x(k); end
                        k = k + 2;
                    otherwise
                        k = k + 1;
                end
            end
            pop(row + j, :) = x;
        end
        row = row + ns;
    end
end

function pop = advisor_seed_init(n_seeds, bounds, var, profile)
% Generate n_seeds from 5 diverse strategy templates built from LLM ratios.
% Each template shifts cruise/high and coast/low1 ratios in a different
% direction so the initial population spans the energy-time Pareto front.
    pop = [];

    if isempty(profile) || ~isstruct(profile) || ~isfield(profile, 'sections')
        return;
    end
    sections = profile.sections;
    if numel(sections) ~= numel(var)
        return;
    end

    dim   = size(bounds, 1);
    lb    = bounds(:,1)';
    ub    = bounds(:,2)';
    rng_b = ub - lb;

    % 5 strategies: [delta_cruise/high, delta_coast/low1]
    %   S1: LLM as-is          S2: speed priority
    %   S3: aggressive energy   S4: high-cruise deep-coast   S5: moderate energy
    deltas = [ 0.00,  0.00;
              +0.10, +0.08;
              -0.10, -0.12;
              +0.06, -0.12;
              -0.05, -0.06 ];
    n_strat = size(deltas, 1);
    n_per   = floor(n_seeds / n_strat);
    counts  = [repmat(n_per, 1, n_strat-1), n_seeds - (n_strat-1)*n_per];

    pop = zeros(n_seeds, dim);
    row = 0;

    for s = 1:n_strat
        dc = deltas(s, 1);
        dr = deltas(s, 2);
        ns = counts(s);

        template = zeros(1, dim);
        k = 1;
        valid = true;
        for i = 1:numel(var)
            sec = sections(i);
            lo1 = lb(k);
            hi1 = ub(k);
            switch var(i)
                case 2
                    if ~isfield(sec,'cruise_ratio') || ~isfield(sec,'coast_ratio')
                        valid = false; break;
                    end
                    cr  = clamp01(double(sec.cruise_ratio) + dc);
                    rr  = clamp01(double(sec.coast_ratio)  + dr);
                    dv1 = lo1 + (hi1 - lo1) * cr;
                    lo2 = lb(k+1); hi2 = min(ub(k+1), dv1);
                    dv2 = lo2 + (hi2 - lo2) * rr;
                    template(k)   = dv1;
                    template(k+1) = min(dv2, dv1);
                    k = k + 2;
                case 3
                    if ~isfield(sec,'high_ratio') || ~isfield(sec,'low1_ratio')
                        valid = false; break;
                    end
                    hr  = clamp01(double(sec.high_ratio)  + dc);
                    lr  = clamp01(double(sec.low1_ratio)  + dr);
                    dv1 = lo1 + (hi1 - lo1) * hr;
                    lo2 = lb(k+1); hi2 = min(ub(k+1), dv1);
                    dv2 = lo2 + (hi2 - lo2) * lr;
                    template(k)   = dv1;
                    template(k+1) = min(dv2, dv1);
                    k = k + 2;
                otherwise
                    template(k) = lo1 + (hi1 - lo1) * clamp01(0.85 + dc);
                    k = k + 1;
            end
        end
        if ~valid, pop = []; return; end

        n_tight = max(1, round(0.6 * ns));
        for j = 1:ns
            noise_scale = 0.025 * (j <= n_tight) + 0.06 * (j > n_tight);
            x = template + (rand(1, dim) - 0.5) .* rng_b * noise_scale;
            x = min(max(x, lb), ub);
            x = enforce_pair_order(x, var);
            pop(row + j, :) = x;
        end
        row = row + ns;
    end
end

function x = enforce_pair_order(x, var)
    k = 1;
    for i = 1:numel(var)
        switch var(i)
            case {2, 3}
                if x(k+1) > x(k)
                    x(k+1) = x(k);
                end
                k = k + 2;
            otherwise
                k = k + 1;
        end
    end
end

function val = clamp01(val)
    val = min(max(double(val), 0.05), 0.99);
end

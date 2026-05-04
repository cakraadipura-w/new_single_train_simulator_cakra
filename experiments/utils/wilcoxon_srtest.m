function p = wilcoxon_srtest(x, y)
%WILCOXON_SRTEST  Wilcoxon signed-rank test (paired).
% Falls back to ranksum if signrank is unavailable (no Statistics Toolbox).
%
% H0: median(x - y) = 0  (paired)

    x = x(:); y = y(:);
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);

    if numel(x) < 2
        p = NaN; return;
    end

    if exist('signrank', 'file') == 2 || exist('signrank', 'builtin') == 5
        p = signrank(x, y);
    elseif exist('ranksum', 'file') == 2
        p = ranksum(x, y);   % unpaired approximation
    else
        % Manual normal approximation of signed-rank statistic
        d = x - y;
        d = d(d ~= 0);
        n = numel(d);
        if n == 0, p = 1; return; end
        [~, idx] = sort(abs(d));
        ranks = zeros(n,1);
        ranks(idx) = 1:n;
        W = sum(ranks(d > 0));
        mu  = n*(n+1)/4;
        sig = sqrt(n*(n+1)*(2*n+1)/24);
        z   = (W - mu) / sig;
        p   = 2 * (1 - normcdf(abs(z)));
    end
end

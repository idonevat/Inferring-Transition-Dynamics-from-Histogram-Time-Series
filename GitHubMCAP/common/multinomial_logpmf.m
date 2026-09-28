function logp = multinomial_logpmf(x, p)
%MULTINOMIAL_LOGPMF Stable multinomial log probability.
x = x(:).'; p = p(:).';
n = sum(x);
if any(p < 0) || abs(sum(p)-1) > 1e-8
    logp = -Inf; return;
end
if any((p==0) & (x>0))
    logp = -Inf; return;
end
idx = x>0;
logp = gammaln(n+1)-sum(gammaln(x+1)) + sum(x(idx).*log(p(idx)));
end

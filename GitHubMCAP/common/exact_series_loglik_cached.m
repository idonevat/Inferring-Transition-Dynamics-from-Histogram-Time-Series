function ll = exact_series_loglik_cached(theta,cache,Delta,K)
%EXACT_SERIES_LOGLIK_CACHED Exact series likelihood with precomputed DP graphs.
[pi_vec,mu]=unpack_theta(theta,K);
P=transition_matrix(pi_vec,mu,Delta);
ll=0;
for t=1:numel(cache)
    z=exact_transition_logprob_cached(cache{t},P);
    if ~isfinite(z), ll=-Inf; return; end
    ll=ll+z;
end
end

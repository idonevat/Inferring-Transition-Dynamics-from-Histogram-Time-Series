function ll = exact_series_loglik(theta,Y,Delta)
%EXACT_SERIES_LOGLIK Exact finite-population observed-data log likelihood.
K=size(Y,2); [pi_vec,mu]=unpack_theta(theta,K);
P=transition_matrix(pi_vec,mu,Delta);
ll=0;
for t=1:size(Y,1)-1
    z=exact_transition_logprob_dp(Y(t,:),Y(t+1,:),P);
    if ~isfinite(z), ll=-Inf; return; end
    ll=ll+z;
end
end

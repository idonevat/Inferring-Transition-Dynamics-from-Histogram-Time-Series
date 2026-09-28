function theta = pack_theta(pi_vec, mu)
%PACK_THETA Convert (pi,mu) to baseline-logit/log unconstrained coordinates.
pi_vec = pi_vec(:).';
mu = mu(:).';
K = numel(pi_vec);
assert(numel(mu) == K-1);
theta_pi = log(pi_vec(1:K-1)./pi_vec(K));
theta_mu = log(mu);
theta = [theta_pi(:); theta_mu(:)];
end

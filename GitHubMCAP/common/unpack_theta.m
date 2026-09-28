function [pi_vec, mu] = unpack_theta(theta, K)
%UNPACK_THETA Baseline-logit coordinates for pi and log coordinates for mu.
theta = theta(:);
assert(numel(theta) == 2*K-2, 'theta must have length 2K-2.');
eta = [theta(1:K-1); 0];
eta = eta - max(eta);
w = exp(eta);
pi_vec = (w/sum(w)).';
mu = exp(theta(K:2*K-2)).';
end

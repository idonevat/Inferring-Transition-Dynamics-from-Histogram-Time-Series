function Q = build_generator(pi_vec, mu)
%BUILD_GENERATOR Reversible nearest-neighbor CTMC generator.
%   Q = BUILD_GENERATOR(pi_vec, mu) constructs the K-by-K generator used in
%   the paper. pi_vec is a positive probability vector and mu has length K-1.

pi_vec = pi_vec(:).';
mu = mu(:).';
K = numel(pi_vec);
assert(numel(mu) == K-1, 'mu must have length K-1.');
assert(all(pi_vec > 0) && abs(sum(pi_vec)-1) < 1e-10, 'pi must be positive and sum to one.');
assert(all(mu > 0), 'mu must be positive.');

Q = zeros(K,K);
for j = 1:K-1
    Q(j,j+1) = mu(j)*sqrt(pi_vec(j+1)/pi_vec(j));
    Q(j+1,j) = mu(j)*sqrt(pi_vec(j)/pi_vec(j+1));
end
Q(1:K+1:end) = -sum(Q,2);
end

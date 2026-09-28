function P = transition_matrix(pi_vec, mu, Delta)
%TRANSITION_MATRIX Finite-time transition matrix exp(Delta Q).
Q = build_generator(pi_vec, mu);
P = expm(Delta*Q);
% Clean very small numerical errors only.
P(abs(P) < 1e-14) = 0;
P = P ./ sum(P,2);
end

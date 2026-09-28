function y_next = simulate_histogram_transition(y, P)
%SIMULATE_HISTOGRAM_TRANSITION Exact aggregate transition under unit CTMC.
y = y(:).';
K = numel(y);
y_next = zeros(1,K);
for j = 1:K
    y_next = y_next + multinomial_sample(y(j), P(j,:));
end
end

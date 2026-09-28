function Y = simulate_histogram_series(y0, P, T)
%SIMULATE_HISTOGRAM_SERIES Simulate T transitions, returning (T+1)-by-K histograms.
y0 = y0(:).';
K = numel(y0);
Y = zeros(T+1,K);
Y(1,:) = y0;
for t = 1:T
    Y(t+1,:) = simulate_histogram_transition(Y(t,:),P);
end
end

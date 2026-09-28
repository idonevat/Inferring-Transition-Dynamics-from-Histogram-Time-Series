function cache = prepare_exact_series_cache(Y)
%PREPARE_EXACT_SERIES_CACHE Precompute exact-DP graphs for all observed transitions.
T=size(Y,1)-1;
cache=cell(T,1);
for t=1:T
    cache{t}=prepare_exact_transition_cache(Y(t,:),Y(t+1,:));
end
end

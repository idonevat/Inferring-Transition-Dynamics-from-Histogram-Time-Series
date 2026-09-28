function logK = exact_transition_logprob_cached(cache,P)
%EXACT_TRANSITION_LOGPROB_CACHED Exact transition probability using cached DP graph.
K=cache.K;
if size(P,1)~=K || size(P,2)~=K, error('P dimension mismatch.'); end

v=1; logscale=0;
for j=1:K
    st=cache.stages{j};
    pj=P(j,:);
    if any(pj<0), logK=-Inf; return; end
    lp=log(max(pj,realmin));
    % Probability for each row composition: multinomial coefficient times p^x.
    logw=st.logcoef + st.X*lp(:);
    w=exp(logw);
    contrib=v(double(st.prev_idx)).*w(double(st.alloc_idx));
    vnext=accumarray(double(st.next_idx),contrib,[st.nnext 1],@sum,0);
    s=sum(vnext);
    if ~(s>0) || ~isfinite(s), logK=-Inf; return; end
    v=vnext/s;
    logscale=logscale+log(s);
end
val=v(cache.target_idx);
if ~(val>0), logK=-Inf; else, logK=log(val)+logscale; end
end

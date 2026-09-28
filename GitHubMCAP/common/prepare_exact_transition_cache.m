function cache = prepare_exact_transition_cache(y,yprime)
%PREPARE_EXACT_TRANSITION_CACHE Precompute parameter-independent exact-DP graph.
% The observed source/destination histograms determine all feasible states,
% row allocations, state transitions, and multinomial coefficients. These are
% cached once and reused for every likelihood evaluation during optimization.

y=round(y(:).'); yprime=round(yprime(:).');
K=numel(y); N=sum(y);
if sum(yprime)~=N || any(y<0) || any(yprime<0)
    error('Source and destination histograms must be nonnegative and have equal totals.');
end
base=N+1;

% states{j+1} contains reachable partial destination histograms after j rows.
states=cell(K+1,1);
states{1}=zeros(1,K);
stages=cell(K,1);
partial_total=0;

for j=1:K
    partial_total=partial_total+y(j);
    X=row_compositions(y(j),K,yprime);
    if isempty(X), error('No feasible row compositions at source row %d.',j); end

    % Multinomial log coefficients n! / prod x_k! are parameter independent.
    logcoef=gammaln(y(j)+1)-sum(gammaln(X+1),2);

    prev=states{j};
    next_map=containers.Map('KeyType','uint64','ValueType','uint32');
    next_states=zeros(0,K);
    prev_idx=zeros(0,1,'uint32');
    alloc_idx=zeros(0,1,'uint32');
    next_idx=zeros(0,1,'uint32');

    for a=1:size(prev,1)
        m=prev(a,:);
        for q=1:size(X,1)
            m2=m+X(q,:);
            if any(m2>yprime), continue; end
            if sum(yprime-m2)~=N-partial_total, continue; end
            key=encode_state(m2,base);
            if isKey(next_map,key)
                ni=next_map(key);
            else
                next_states(end+1,:)=m2; %#ok<AGROW>
                ni=uint32(size(next_states,1));
                next_map(key)=ni;
            end
            prev_idx(end+1,1)=uint32(a); %#ok<AGROW>
            alloc_idx(end+1,1)=uint32(q); %#ok<AGROW>
            next_idx(end+1,1)=ni; %#ok<AGROW>
        end
    end

    stages{j}=struct('X',X,'logcoef',logcoef,'prev_idx',prev_idx, ...
        'alloc_idx',alloc_idx,'next_idx',next_idx,'nnext',size(next_states,1));
    states{j+1}=next_states;
end

% Final feasible state should be yprime; retain its location defensively.
idx=find(all(states{K+1}==yprime,2),1);
if isempty(idx), error('Target histogram is not reachable in cached DP graph.'); end
cache=struct('y',y,'yprime',yprime,'K',K,'N',N,'stages',{stages}, ...
    'target_idx',idx,'states',{states});
end

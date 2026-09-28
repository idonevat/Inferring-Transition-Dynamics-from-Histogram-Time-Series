function logK = exact_transition_logprob_dp(y, yprime, P)
%EXACT_TRANSITION_LOGPROB_DP Exact aggregate transition log probability.
% Implements the forward recursion in Appendix app:exactcomp.

y = round(y(:).'); yprime = round(yprime(:).');
K = numel(y); N = sum(y);
if sum(yprime) ~= N || any(y<0) || any(yprime<0)
    logK = -Inf; return;
end
if size(P,1)~=K || size(P,2)~=K
    error('P dimension mismatch.');
end
base = N+1;
% Map state key -> probability. Scaling after each row avoids underflow.
A = containers.Map('KeyType','uint64','ValueType','double');
A(encode_state(zeros(1,K),base)) = 1;
logscale = 0;
partial_total = 0;
for j = 1:K
    partial_total = partial_total + y(j);
    X = row_compositions(y(j),K,yprime);
    rowprob = zeros(size(X,1),1);
    for q = 1:size(X,1)
        lp = multinomial_logpmf(X(q,:),P(j,:));
        rowprob(q) = exp(lp);
    end
    B = containers.Map('KeyType','uint64','ValueType','double');
    keysA = A.keys; valsA = A.values;
    for ia = 1:numel(keysA)
        keyA = keysA{ia}; valA = valsA{ia};
        % Decode is avoided: store decoded vector via base representation.
        m = decode_state_local(keyA,base,K);
        for q = 1:size(X,1)
            if rowprob(q)==0, continue; end
            m2 = m + X(q,:);
            if any(m2 > yprime), continue; end
            % Feasibility: remaining source mass must be able to fill target.
            if sum(yprime-m2) ~= N-partial_total, continue; end
            key2 = encode_state(m2,base);
            add = valA*rowprob(q);
            if isKey(B,key2), B(key2)=B(key2)+add; else, B(key2)=add; end
        end
    end
    if B.Count==0, logK=-Inf; return; end
    vb = cell2mat(B.values);
    s = sum(vb);
    kb = B.keys;
    for ii=1:numel(kb), B(kb{ii}) = B(kb{ii})/s; end
    logscale = logscale + log(s);
    A = B;
end
ktarget = encode_state(yprime,base);
if ~isKey(A,ktarget), logK=-Inf; else, logK=log(A(ktarget))+logscale; end
end

function x = decode_state_local(key,base,K)
x = zeros(1,K); q = uint64(key); b=uint64(base);
for k=1:K
    x(k)=double(mod(q,b));
    q=idivide(q,b,'floor');
end
end

function X = row_compositions(n, K, upper)
%ROW_COMPOSITIONS All K-part compositions of n bounded componentwise by upper.
% Intended for exact DP at small/moderate N and K.
if nargin < 3 || isempty(upper), upper = n*ones(1,K); end
upper = floor(upper(:).');
rows = {};
cur = zeros(1,K);
recurse(1,n);
X = vertcat(rows{:});
if isempty(X), X = zeros(0,K); end

    function recurse(k, remaining)
        if k == K
            if remaining <= upper(k)
                cur(k) = remaining;
                rows{end+1,1} = cur; %#ok<AGROW>
            end
            return;
        end
        maxx = min(remaining,upper(k));
        for v = 0:maxx
            cur(k)=v;
            recurse(k+1,remaining-v);
        end
    end
end

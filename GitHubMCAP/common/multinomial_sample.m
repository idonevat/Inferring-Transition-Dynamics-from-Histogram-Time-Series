function counts = multinomial_sample(n, p)
%MULTINOMIAL_SAMPLE Toolbox-free multinomial draw.
p = p(:).';
p = max(p,0); p = p/sum(p);
K = numel(p);
counts = zeros(1,K);
if n == 0, return; end
u = rand(n,1);
c = cumsum(p); c(end)=1;
for i = 1:n
    k = find(u(i) <= c,1,'first');
    counts(k) = counts(k)+1;
end
end

function y = make_source_histogram(N,rho)
%MAKE_SOURCE_HISTOGRAM Deterministic integer histogram close to N*rho.
rho=max(rho(:).',0); rho=rho/sum(rho);
z=N*rho; y=floor(z); rem=N-sum(y);
[~,ord]=sort(z-y,'descend');
y(ord(1:rem))=y(ord(1:rem))+1;
end

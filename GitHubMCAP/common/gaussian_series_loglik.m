function ll = gaussian_series_loglik(theta,Y,Delta)
%GAUSSIAN_SERIES_LOGLIK Large-population Gaussian log likelihood (constants omitted).
K=size(Y,2); d=K-1; [pi_vec,mu]=unpack_theta(theta,K);
P=transition_matrix(pi_vec,mu,Delta);
N=sum(Y(1,:)); ll=0;
for t=1:size(Y,1)-1
    rho=Y(t,:)/N;
    q=P(:,1:d); % row j reduced destination vector
    m=(rho*P); m=m(1:d);
    V=zeros(d,d);
    for j=1:K
        qj=q(j,:).';
        V=V+rho(j)*(diag(qj)-qj*qj.');
    end
    V=(V+V.')/2;
    % Numerical regularization only if necessary.
    [R,pflag]=chol(V);
    if pflag~=0
        jitter=max(1e-12,1e-10*trace(V)/max(d,1));
        V=V+jitter*eye(d); [R,pflag]=chol(V);
    end
    if pflag~=0, ll=-Inf; return; end
    r=Y(t+1,1:d)/N-m;
    logdet=2*sum(log(diag(R)));
    quad=(r/R)*(r/R).'; % r V^{-1} r'
    ll=ll-0.5*(logdet+N*quad);
end
end

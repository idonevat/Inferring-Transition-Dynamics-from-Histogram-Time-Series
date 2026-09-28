function fit = fit_model(Y,Delta,method,start_theta,nstarts,varargin)
%FIT_MODEL Fit exact or Gaussian likelihood using fminsearch.
% Optional sixth argument: a precomputed exact-series cache for method
% 'exact_cached'. Passing it avoids rebuilding the combinatorial cache.
if nargin<5 || isempty(nstarts), nstarts=3; end
K=size(Y,2);
if nargin<4 || isempty(start_theta)
    p0=mean(Y,1); p0=max(p0,1); p0=p0/sum(p0);
    m0=0.2*ones(1,K-1);
    start_theta=pack_theta(p0,m0);
end
switch lower(method)
    case 'exact', loglik=@(th) exact_series_loglik(th,Y,Delta);
    case 'exact_cached'
        if ~isempty(varargin) && ~isempty(varargin{1})
            exact_cache=varargin{1};
        else
            exact_cache=prepare_exact_series_cache(Y);
        end
        loglik=@(th) exact_series_loglik_cached(th,exact_cache,Delta,K);
    case 'gaussian', loglik=@(th) gaussian_series_loglik(th,Y,Delta);
    otherwise, error('method must be exact, exact_cached, or gaussian');
end
bestf=Inf; bestth=[]; bestexit=[]; bestoutput=[];
opts=optimset('Display','off','MaxFunEvals',5000,'MaxIter',2000,'TolX',1e-6,'TolFun',1e-6);
for s=1:nstarts
    if s==1, th0=start_theta(:); else, th0=start_theta(:)+0.25*randn(size(start_theta(:))); end
    obj=@(th) safe_negloglik(th,loglik);
    [th,fval,exitflag,output]=fminsearch(obj,th0,opts);
    if fval<bestf
        bestf=fval; bestth=th; bestexit=exitflag; bestoutput=output;
    end
end
[pi_hat,mu_hat]=unpack_theta(bestth,K);
fit.theta=bestth; fit.pi=pi_hat; fit.mu=mu_hat;
fit.P=transition_matrix(pi_hat,mu_hat,Delta);
fit.loglik=-bestf; fit.exitflag=bestexit; fit.output=bestoutput;
end

% function f=safe_negloglik(th,loglik)
% ll=loglik(th);
% if ~isfinite(ll), f=1e100; else, f=-ll; end
% end
function nll = safe_negloglik(th,loglik)

try
    ll = loglik(th);

    if isempty(ll) || ~isscalar(ll) || ~isreal(ll) || ~isfinite(ll)
        nll = Inf;
    else
        nll = -ll;
    end

catch
    % fminsearch may explore numerically invalid parameter values
    % (e.g. exp-transformed mobility underflowing to zero).
    % Treat such points as outside the admissible parameter space.
    nll = Inf;
end

end
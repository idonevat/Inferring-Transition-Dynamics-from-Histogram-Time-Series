% Public reproducibility script.
% Expected final result file: experiment4_stage1_population_results_R50.mat
% After running, move/copy the result into the matching results/experiment*/ folder
% if the script saves it in the current MATLAB working directory.

%% experiment4_stage1_population_full.m
% EXPERIMENT 4A (FINAL): Population density and state-space scalability.
%
% Design: K={4,10,20,50,100}, N/K={10,50,200,1000}, T=100, R=50.
% Exact aggregate simulation; Gaussian aggregate likelihood fitting.
% Every fit is independent and uses mobility starts 0.10, 0.25, 0.50.

clear; clc;
K_grid=[4 10 20 50 100];
occ_grid=[10 50 200 1000];
T=100; Delta=1; R=50;
mu_start_levels=[0.10 0.25 0.50];
base_seed=20260907;
max_iterations=700; max_function_evals=4500;
do_gradient_check=true;
checkpoint_file='experiment4_stage1_population_checkpoint_R50.mat';
results_file='experiment4_stage1_population_results_R50.mat';

if exist('fminunc','file')~=2, error('This experiment requires fminunc.'); end
fprintf('\nExperiment 4A final: population density and state-space scalability\n');
fprintf('R=%d, T=%d, K=%s, N/K=%s\n\n',R,T,mat2str(K_grid),mat2str(occ_grid));

if do_gradient_check
    rng(base_seed,'twister'); Kc=6; Nc=300; Tc=8;
    [pic,muc,rhoc]=make_profiles(Kc); Pc=transition_matrix_local(pic,muc,Delta);
    Yc=simulate_histogram_series_fast(make_source_histogram_local(Nc,rhoc),Pc,Tc);
    th=pack_theta_local(pic,muc)+0.03*randn(2*Kc-2,1);
    check_gaussian_gradient(th,Yc,Delta); fprintf('Gradient check passed.\n\n');
end

nK=numel(K_grid); nO=numel(occ_grid);
if exist(checkpoint_file,'file')
    S=load(checkpoint_file);
    if ~isequal(S.K_grid,K_grid)||~isequal(S.occ_grid,occ_grid)||S.T~=T||S.Delta~=Delta||S.R~=R||~isequal(S.mu_start_levels,mu_start_levels)
        error('Existing Stage-1 checkpoint has a different design.');
    end
    raw=S.raw; completed=S.completed;
    fprintf('Resuming checkpoint: %d/%d fits complete.\n\n',nnz(completed),numel(completed));
else
    raw.err_pi=nan(nK,nO,R); raw.err_mu=nan(nK,nO,R); raw.err_Q=nan(nK,nO,R); raw.err_P=nan(nK,nO,R);
    raw.fit_time=nan(nK,nO,R); raw.exitflag=nan(nK,nO,R); raw.iterations=nan(nK,nO,R);
    raw.final_nll=nan(nK,nO,R); raw.start_index=nan(nK,nO,R); raw.max_mu_hat=nan(nK,nO,R);
    raw.pi_hat=cell(nK,nO); raw.mu_hat=cell(nK,nO); raw.true_pi=cell(nK,1); raw.true_mu=cell(nK,1); raw.true_rho=cell(nK,1);
    for iK=1:nK, K=K_grid(iK); for iO=1:nO, raw.pi_hat{iK,iO}=nan(R,K); raw.mu_hat{iK,iO}=nan(R,K-1); end, end
    completed=false(nK,nO,R);
end

for iK=1:nK
    K=K_grid(iK); [pi_true,mu_true,rho_true]=make_profiles(K);
    P_true=transition_matrix_local(pi_true,mu_true,Delta); Q_true=build_generator_local(pi_true,mu_true);
    raw.true_pi{iK}=pi_true; raw.true_mu{iK}=mu_true; raw.true_rho{iK}=rho_true;
    for iO=1:nO
        occ=occ_grid(iO); N=K*occ;
        for r=1:R
            if completed(iK,iO,r), continue; end
            rng(base_seed+100000*iK+1000*iO+r,'twister');
            y0=make_source_histogram_local(N,rho_true); Y=simulate_histogram_series_fast(y0,P_true,T);
            p0=mean(Y,1)+0.5; p0=p0/sum(p0); mu_starts=repmat(mu_start_levels(:),1,K-1);
            tic; fit=fit_gaussian_multistart(Y,Delta,p0,mu_starts,max_iterations,max_function_evals); raw.fit_time(iK,iO,r)=toc;
            pi_hat=fit.pi; mu_hat=fit.mu; P_hat=fit.P; Q_hat=build_generator_local(pi_hat,mu_hat);
            raw.err_pi(iK,iO,r)=norm(pi_hat-pi_true)/norm(pi_true); raw.err_mu(iK,iO,r)=norm(mu_hat-mu_true)/norm(mu_true);
            raw.err_Q(iK,iO,r)=norm(Q_hat-Q_true,'fro')/norm(Q_true,'fro'); raw.err_P(iK,iO,r)=norm(P_hat-P_true,'fro')/norm(P_true,'fro');
            raw.pi_hat{iK,iO}(r,:)=pi_hat; raw.mu_hat{iK,iO}(r,:)=mu_hat; raw.max_mu_hat(iK,iO,r)=max(mu_hat);
            raw.exitflag(iK,iO,r)=fit.exitflag; raw.iterations(iK,iO,r)=fit.output.iterations; raw.final_nll(iK,iO,r)=-fit.loglik; raw.start_index(iK,iO,r)=fit.start_index;
            completed(iK,iO,r)=true;
            save(checkpoint_file,'K_grid','occ_grid','T','Delta','R','mu_start_levels','raw','completed','base_seed','max_iterations','max_function_evals','-v7.3');
            fprintf('K=%3d N/K=%4d N=%7d rep=%2d/%2d | E_mu=%7.3g E_P=%7.3g maxmu=%8.3g time=%7.1fs exit=%2d start=%d\n',K,occ,N,r,R,raw.err_mu(iK,iO,r),raw.err_P(iK,iO,r),raw.max_mu_hat(iK,iO,r),raw.fit_time(iK,iO,r),fit.exitflag,fit.start_index);
        end
    end
end

results.K_grid=K_grid; results.occ_grid=occ_grid; results.T=T; results.Delta=Delta; results.R=R; results.mu_start_levels=mu_start_levels; results.raw=raw;
results.med_err_pi=median(raw.err_pi,3,'omitnan'); results.med_err_mu=median(raw.err_mu,3,'omitnan'); results.med_err_Q=median(raw.err_Q,3,'omitnan'); results.med_err_P=median(raw.err_P,3,'omitnan');
results.q90_err_mu=prctile(raw.err_mu,90,3); results.med_time=median(raw.fit_time,3,'omitnan'); results.exit1_rate=mean(raw.exitflag==1,3,'omitnan'); results.positive_exit_rate=mean(raw.exitflag>0,3,'omitnan'); results.runaway10_rate=mean(raw.max_mu_hat>10,3,'omitnan');
save(results_file,'results','-v7.3');
fprintf('\nMedian relative mu error:\n'); disp(array2table(results.med_err_mu,'VariableNames',compose('occ_%d',occ_grid),'RowNames',compose('K_%d',K_grid)));
fprintf('90th percentile relative mu error:\n'); disp(array2table(results.q90_err_mu,'VariableNames',compose('occ_%d',occ_grid),'RowNames',compose('K_%d',K_grid)));
fprintf('Median fitting time (seconds):\n'); disp(array2table(results.med_time,'VariableNames',compose('occ_%d',occ_grid),'RowNames',compose('K_%d',K_grid)));
fprintf('Fraction with max(mu_hat)>10:\n'); disp(array2table(results.runaway10_rate,'VariableNames',compose('occ_%d',occ_grid),'RowNames',compose('K_%d',K_grid)));
fprintf('Fraction with exitflag==1:\n'); disp(array2table(results.exit1_rate,'VariableNames',compose('occ_%d',occ_grid),'RowNames',compose('K_%d',K_grid)));
fprintf('\nSaved %s\n',results_file);

%% LOCAL FUNCTIONS
function [pi_vec,mu,rho] = make_profiles(K)
% Smooth profiles defined on a common ordered state domain.
x  = ((1:K)-0.5)/K;
xb = (1:K-1)/K;

wpi = 1 + 0.35*cos(2*pi*x) + 0.20*sin(4*pi*x);
if any(wpi<=0), error('Equilibrium profile is not positive.'); end
pi_vec = wpi/sum(wpi);

mu = 0.25 + 0.08*sin(2*pi*xb) + 0.04*cos(4*pi*xb);
if any(mu<=0), error('Mobility profile is not positive.'); end

% Moderate, smooth departure from equilibrium: enough mean excitation without
% making initialization unrealistically difficult.
wrho = pi_vec .* exp(0.30*sin(2*pi*x) - 0.18*cos(4*pi*x));
rho = wrho/sum(wrho);
end

function Q = build_generator_local(pi_vec,mu)
pi_vec = pi_vec(:).';
mu = mu(:).';
K = numel(pi_vec);
Q = zeros(K,K);
for j=1:K-1
    a = mu(j)*sqrt(pi_vec(j+1)/pi_vec(j));
    b = mu(j)*sqrt(pi_vec(j)/pi_vec(j+1));
    Q(j,j+1)=a;
    Q(j+1,j)=b;
end
Q(1:K+1:end) = -sum(Q,2);
end

function P = transition_matrix_local(pi_vec,mu,Delta)
P = expm(Delta*build_generator_local(pi_vec,mu));
P(P<0 & P>-1e-13)=0;
P = P./sum(P,2);
end

function theta = pack_theta_local(pi_vec,mu)
K = numel(pi_vec);
theta = [log(pi_vec(1:K-1)./pi_vec(K)), log(mu(:).')].';
end

function [pi_vec,mu] = unpack_theta_local(theta,K)
theta = theta(:);
if numel(theta)~=2*K-2
    error('theta has wrong dimension.');
end
eta = [theta(1:K-1);0];
eta = eta-max(eta);
w = exp(eta);
pi_vec = (w/sum(w)).';
mu = exp(theta(K:2*K-2)).';
end

function y = make_source_histogram_local(N,p)
% Largest-remainder deterministic rounding, preserving total N.
a = N*p(:).';
y = floor(a);
left = N-sum(y);
[~,ord] = sort(a-y,'descend');
if left>0
    y(ord(1:left)) = y(ord(1:left))+1;
end
end

function Y = simulate_histogram_series_fast(y0,P,T)
y0 = y0(:).';
K = numel(y0);
Y = zeros(T+1,K);
Y(1,:) = y0;
for t=1:T
    yn = zeros(1,K);
    for j=1:K
        if Y(t,j)>0
            yn = yn + multinomial_sample_fast(Y(t,j),P(j,:));
        end
    end
    Y(t+1,:) = yn;
end
end

function counts = multinomial_sample_fast(n,p)
% Exact multinomial draw. Sequential binomials are O(K) when binornd is
% available; otherwise a vectorized inverse-CDF draw is used.
p = max(p(:).',0);
p = p/sum(p);
K = numel(p);
counts = zeros(1,K);
if n==0, return; end

if exist('binornd','file')==2
    remain_n = n;
    remain_p = 1;
    for k=1:K-1
        if remain_n==0, break; end
        pk = min(max(p(k)/remain_p,0),1);
        counts(k) = binornd(remain_n,pk);
        remain_n = remain_n-counts(k);
        remain_p = remain_p-p(k);
        if remain_p<=eps, break; end
    end
    counts(K)=remain_n;
else
    edges = [0 cumsum(p)];
    edges(end)=1;
    counts = histcounts(rand(n,1),edges);
end
end

function fit = fit_gaussian_multistart(Y,Delta,pi0,mu_starts,maxit,maxfev)
% Independent explicit multistart optimization for one fixed dataset.
% Each row of mu_starts is either scalar (expanded across interfaces) or a
% complete K-1 vector. The best finite Gaussian negative log-likelihood wins.

K=size(Y,2);

opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'SpecifyObjectiveGradient',true, ...
    'Display','off', ...
    'MaxIterations',maxit, ...
    'MaxFunctionEvaluations',maxfev, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-8, ...
    'FunctionTolerance',1e-8);

bestf=Inf; best=[];
nstarts=size(mu_starts,1);

for s=1:nstarts
    ms=mu_starts(s,:);
    if isscalar(ms)
        m0=ms*ones(1,K-1);
    elseif numel(ms)==K-1
        m0=reshape(ms,1,K-1);
    else
        error('mu_starts row has incompatible dimension.');
    end

    th0=pack_theta_local(pi0,m0);
    fun=@(th) gaussian_nll_score(th,Y,Delta);
    [th,fval,exitflag,output]=fminunc(fun,th0,opts);

    if isfinite(fval) && fval<bestf
        bestf=fval;
        best.theta=th;
        best.exitflag=exitflag;
        best.output=output;
        best.start_index=s;
        best.start_mu=m0;
    end
end

if isempty(best)
    error('All Gaussian optimization starts failed.');
end

[pi_hat,mu_hat]=unpack_theta_local(best.theta,K);
fit.theta=best.theta;
fit.pi=pi_hat;
fit.mu=mu_hat;
fit.P=transition_matrix_local(pi_hat,mu_hat,Delta);
fit.loglik=-bestf;
fit.exitflag=best.exitflag;
fit.output=best.output;
fit.start_index=best.start_index;
fit.start_mu=best.start_mu;
end

function [nll,g] = gaussian_nll_score(theta,Y,Delta)
% Gaussian aggregate negative log likelihood and analytic gradient.
%
% The likelihood is the same reduced-coordinate Gaussian likelihood used in
% gaussian_series_loglik.m. Constants independent of theta are omitted.

K = size(Y,2);
d = K-1;
N = sum(Y(1,:));

[pi_vec,mu] = unpack_theta_local(theta,K);
if any(~isfinite(pi_vec)) || any(~isfinite(mu)) || any(mu<=0) || max(mu)>1e6
    nll=1e100;
    g=zeros(size(theta));
    return;
end

Q = build_generator_local(pi_vec,mu);
A = Delta*Q;
P = expm(A);

if any(~isfinite(P(:)))
    nll=1e100;
    g=zeros(size(theta));
    return;
end

q = P(:,1:d);
ll = 0;
GP = zeros(K,K);

for t=1:size(Y,1)-1
    rho = Y(t,:)/N;
    mfull = rho*P;
    m = mfull(1:d).';

    % V = sum_j rho_j[diag(q_j)-q_j q_j']
    V = diag(m) - q.'*(rho(:).*q);
    V = (V+V.')/2;

    [R,pflag] = chol(V);  % V = R'*R
    if pflag~=0
        jitter=max(1e-12,1e-10*trace(V)/max(d,1));
        V=V+jitter*eye(d);
        [R,pflag]=chol(V);
    end
    if pflag~=0
        nll=1e100;
        g=zeros(size(theta));
        return;
    end

    r = Y(t+1,1:d).'/N - m;
    z = R\(R'\r);
    logdet = 2*sum(log(diag(R)));
    quad = r.'*z;
    ll = ll - 0.5*(logdet + N*quad);

    % Differential:
    % dl = (N z)' dm + tr(GV' dV),
    % GV = .5*(N z z' - V^{-1}).
    Vinv = R\(R'\eye(d));
    GV = 0.5*(N*(z*z.') - Vinv);
    GV = (GV+GV.')/2;
    gm = N*z;
    dg = diag(GV);

    % Gradient with respect to the first d entries of each row of P.
    % Each row j is weighted by source proportion rho_j.
    common = ones(K,1)*(gm+dg).' - 2*q*GV;
    GP(:,1:d) = GP(:,1:d) + rho(:).*common;
end

% Adjoint Frechet derivative through P = expm(A).
% If dP = L_A(dA), then grad_A = L_{A'}(grad_P).
Z = zeros(K,K);
B = [A.', GP; Z, A.'];
EB = expm(B);
GA = EB(1:K,K+1:2*K);
GQ = Delta*GA;

% Chain Q -> log(pi), log(mu).
g_s = zeros(K,1);
g_logmu = zeros(K-1,1);

for j=1:K-1
    a = Q(j,j+1);
    b = Q(j+1,j);

    % derivative with respect to log(mu_j)
    g_logmu(j) = ...
        a*(GQ(j,j+1)-GQ(j,j)) + ...
        b*(GQ(j+1,j)-GQ(j+1,j+1));

    % derivative wrt s_j=log(pi_j) and s_{j+1}=log(pi_{j+1})
    c = 0.5*a*(GQ(j,j)-GQ(j,j+1)) + ...
        0.5*b*(GQ(j+1,j)-GQ(j+1,j+1));
    g_s(j)   = g_s(j)   + c;
    g_s(j+1) = g_s(j+1) - c;
end

% theta_pi are baseline logits eta_1,...,eta_{K-1}, eta_K=0.
g_eta = g_s(1:K-1) - pi_vec(1:K-1).'*sum(g_s);
g_ll = [g_eta; g_logmu];

nll = -ll;
g = -g_ll;

if ~isfinite(nll) || any(~isfinite(g))
    nll=1e100;
    g=zeros(size(theta));
end
end

function check_gaussian_gradient(theta,Y,Delta)
[f,g] = gaussian_nll_score(theta,Y,Delta);
if ~isfinite(f), error('Gradient-check point has nonfinite objective.'); end

rng(9817,'twister');
m = numel(theta);
idx = unique(round(linspace(1,m,min(m,10))));
idx = unique([idx, randi(m,1,min(5,m))]);

rel = zeros(size(idx));
for ii=1:numel(idx)
    j=idx(ii);
    h=1e-6*(1+abs(theta(j)));
    ep=zeros(m,1); ep(j)=h;
    fp=gaussian_nll_score(theta+ep,Y,Delta);
    fm=gaussian_nll_score(theta-ep,Y,Delta);
    gn=(fp-fm)/(2*h);
    rel(ii)=abs(gn-g(j))/max([1,abs(gn),abs(g(j))]);
end
fprintf('  max relative coordinate-gradient discrepancy = %.3e\n',max(rel));
if max(rel)>5e-4
    error('Analytic Gaussian gradient failed finite-difference check.');
end
end

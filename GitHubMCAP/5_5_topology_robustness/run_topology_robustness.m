%% experiment_nonlocal_tree_robustness.m
% Robustness check: path graph versus a nonlocal branching tree.
%
% Both models have K-1 edges and therefore the SAME number of mobility
% parameters.  The tree is not a nearest-neighbour path.  The purpose is not
% to declare one graph easier, but to check whether the qualitative inference
% behavior seen in the paper persists beyond nearest-neighbour movement.
%
% Data are generated exactly from the unit-level CTMC and only aggregate
% histograms are fitted, using the same Gaussian aggregate likelihood and
% analytic-gradient strategy as Experiment 4.
%
% Requires Optimization Toolbox (fminunc).

clear; clc;

%% ---------------- USER SETTINGS -----------------------------------------
RUN_MODE = 'full';              % 'pilot' or 'full'
K = 10;
N = 100;                         % N/K = 10
Delta = 1;
T_grid = [20 50 100 200];
Tmax = max(T_grid);

switch lower(RUN_MODE)
    case 'pilot'
        Rmc = 8;
    case 'full'
        Rmc = 50;
    otherwise
        error('RUN_MODE must be ''pilot'' or ''full''.');
end

base_seed = 20260928;
nstarts = 2;
max_iterations = 600;
max_function_evals = 4000;
do_gradient_check = true;

results_file = sprintf('nonlocal_tree_robustness_%s.mat',lower(RUN_MODE));

if exist('fminunc','file') ~= 2
    error('This experiment requires fminunc (Optimization Toolbox).');
end

%% ---------------- GRAPH DEFINITIONS -------------------------------------
% Path: the numerical specialization used in the paper.
E_path = [(1:K-1)' (2:K)'];

% Branching tree: same K-1 edges/parameters, but not a path and includes
% direct transitions between non-consecutive state labels.
E_tree = [1 2; 1 3; 2 4; 2 5; 3 6; 3 7; 4 8; 5 9; 7 10];

assert(size(E_path,1)==K-1 && size(E_tree,1)==K-1);
assert(is_connected_tree(E_tree,K),'E_tree must be a connected tree.');

graph_names = {'path','tree'};
graph_edges = {E_path,E_tree};
G = numel(graph_names);

%% ---------------- COMMON TRUE PARAMETERS --------------------------------
% Smooth positive equilibrium profile and moderate non-equilibrium start.
x = ((1:K)-0.5)/K;
wpi = 1 + 0.35*cos(2*pi*x) + 0.20*sin(4*pi*x);
pi_true = wpi/sum(wpi);

% Same 9 edge-mobility values are used in both graphs.  Only graph topology
% changes.  Values are in the same range as the paper's numerical studies.
e = (1:K-1)/(K-1);
mu_true = 0.25 + 0.07*sin(2*pi*e) + 0.035*cos(4*pi*e);

wrho = pi_true .* exp(0.30*sin(2*pi*x) - 0.18*cos(4*pi*x));
rho0 = wrho/sum(wrho);
y0 = make_source_histogram_local(N,rho0);

fprintf('\nNon-nearest-neighbour robustness experiment\n');
fprintf('Mode=%s, R=%d, K=%d, N=%d, Delta=%g\n',RUN_MODE,Rmc,K,N,Delta);
fprintf('T = %s\n',mat2str(T_grid));
fprintf('Both graphs have %d mobility parameters.\n\n',K-1);

%% ---------------- STORAGE ------------------------------------------------
nT = numel(T_grid);
err_pi = nan(G,nT,Rmc);
err_mu = nan(G,nT,Rmc);
err_P  = nan(G,nT,Rmc);
fit_time = nan(G,nT,Rmc);
max_mu = nan(G,nT,Rmc);
exitflag = nan(G,nT,Rmc);

%% ---------------- GRADIENT CHECK ----------------------------------------
if do_gradient_check
    for gg=1:G
        E = graph_edges{gg};
        Ptrue = transition_matrix_graph(pi_true,mu_true,E,Delta);
        rng(base_seed + 1000*gg,'twister');
        Ytest = simulate_histogram_series_fast(y0,Ptrue,30);
        thtrue = pack_theta_graph(pi_true,mu_true);
        fprintf('Gradient check: %s graph\n',graph_names{gg});
        check_gaussian_gradient_graph(thtrue,Ytest,E,Delta);
    end
end

%% ---------------- MONTE CARLO -------------------------------------------
for gg=1:G
    E = graph_edges{gg};
    Qtrue = build_generator_graph(pi_true,mu_true,E);
    Ptrue = expm(Delta*Qtrue);
    Ptrue(Ptrue<0 & Ptrue>-1e-13)=0;
    Ptrue = Ptrue./sum(Ptrue,2);

    fprintf('\n--- %s graph ---\n',upper(graph_names{gg}));

    for rr=1:Rmc
        % One Tmax trajectory per replication; shorter T values are prefixes.
        rng(base_seed + 100000*gg + rr,'twister');
        Yfull = simulate_histogram_series_fast(y0,Ptrue,Tmax);

        for it=1:nT
            T = T_grid(it);
            Y = Yfull(1:T+1,:);

            % Start near truth, but not at truth.  Same perturbation scale for
            % both graph families.  Additional starts are generated internally.
            rng(base_seed + 1000000*gg + 1000*rr + it,'twister');
            pi0 = normalize_positive(pi_true .* exp(0.10*randn(1,K)));
            mu0 = mu_true .* exp(0.15*randn(1,K-1));
            th0 = pack_theta_graph(pi0,mu0);

            tic;
            fit = fit_gaussian_graph(Y,E,Delta,th0,nstarts, ...
                max_iterations,max_function_evals);
            fit_time(gg,it,rr) = toc;

            err_pi(gg,it,rr) = norm(fit.pi-pi_true,2);
            err_mu(gg,it,rr) = norm(fit.mu-mu_true,2)/norm(mu_true,2);
            err_P(gg,it,rr)  = norm(fit.P-Ptrue,'fro')/norm(Ptrue,'fro');
            max_mu(gg,it,rr) = max(fit.mu);
            exitflag(gg,it,rr) = fit.exitflag;
        end

        fprintf('  replication %d/%d complete\n',rr,Rmc);
    end
end

%% ---------------- SUMMARIES ---------------------------------------------
med_pi = median(err_pi,3,'omitnan');
med_mu = median(err_mu,3,'omitnan');
med_P  = median(err_P,3,'omitnan');
med_time = median(fit_time,3,'omitnan');
pct_large = 100*mean(max_mu>10,3,'omitnan');

fprintf('\nMedian relative mobility error\n');
for gg=1:G
    fprintf('  %-5s: %s\n',graph_names{gg},mat2str(med_mu(gg,:),4));
end
fprintf('Median relative P error\n');
for gg=1:G
    fprintf('  %-5s: %s\n',graph_names{gg},mat2str(med_P(gg,:),4));
end
fprintf('Percent fits with max(mu_hat)>10\n');
for gg=1:G
    fprintf('  %-5s: %s\n',graph_names{gg},mat2str(pct_large(gg,:),4));
end
fprintf('Median fit time (s)\n');
for gg=1:G
    fprintf('  %-5s: %s\n',graph_names{gg},mat2str(med_time(gg,:),4));
end

%% ---------------- SAVE ---------------------------------------------------
results = struct;
results.RUN_MODE = RUN_MODE;
results.K = K; results.N = N; results.Delta = Delta;
results.T_grid = T_grid; results.R = Rmc;
results.pi_true = pi_true; results.mu_true = mu_true; results.rho0 = rho0;
results.E_path = E_path; results.E_tree = E_tree;
results.graph_names = graph_names;
results.err_pi = err_pi; results.err_mu = err_mu; results.err_P = err_P;
results.fit_time = fit_time; results.max_mu = max_mu; results.exitflag = exitflag;
results.med_pi = med_pi; results.med_mu = med_mu; results.med_P = med_P;
results.med_time = med_time; results.pct_large = pct_large;
save(results_file,'results','-v7.3');

%% plot_nonlocal_tree_robustness.m
% Plot the full nonlocal-tree robustness experiment.
%
% Expected input:
%   nonlocal_tree_robustness_full.mat
%
% Figure:
%   (a) equilibrium composition error E_pi
%   (b) mobility relative error E_mu
%   (c) transition-matrix error E_P
%
% Each panel shows the median across Monte Carlo replications, with
% 25th--75th percentile (IQR) error bars, for the nearest-neighbour and nonlocal models.
%
% The plotting order matches the numerical figures in the paper:
% E_pi -> E_mu -> E_P.

clear; clc; close all;

%% ---------------- SETTINGS ----------------------------------------------
results_file = 'nonlocal_tree_robustness_full.mat';

% Output names
fig_file = 'nonlocal_tree_robustness_figure.fig';
pdf_file = 'nonlocal_tree_robustness_figure.pdf';
png_file = 'nonlocal_tree_robustness_figure.png';

%% ---------------- LOAD RESULTS ------------------------------------------
S = load(results_file);

if isfield(S,'results')
    R = S.results;
else
    R = S;
end

T_grid = R.T_grid(:);

% Dimensions in saved file:
%   err_*(replication, T index, graph)
% graph 1 = nearest-neighbour
% graph 2 = nonlocal
err_pi = R.err_pi;
err_mu = R.err_mu;
err_P  = R.err_P;

graph_labels = {'Nearest-neighbour','Nonlocal'};

%% ---------------- SUMMARIES ---------------------------------------------
% median and IQR across replications (dimension 1)
med_pi = squeeze(median(err_pi,1,'omitnan'));
med_mu = squeeze(median(err_mu,1,'omitnan'));
med_P  = squeeze(median(err_P ,1,'omitnan'));

q25_pi = squeeze(prctile(err_pi,25,1));
q75_pi = squeeze(prctile(err_pi,75,1));

q25_mu = squeeze(prctile(err_mu,25,1));
q75_mu = squeeze(prctile(err_mu,75,1));

q25_P = squeeze(prctile(err_P,25,1));
q75_P = squeeze(prctile(err_P,75,1));

%% ---------------- FIGURE ------------------------------------------------
fig = figure('Color','w','Position',[100 100 1180 350]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

lw = 1.6;
ms = 6;

% MATLAB default axes color order; use the first two entries consistently.
co = get(groot,'defaultAxesColorOrder');
c1 = co(1,:);
c2 = co(2,:);

% Helper for an IQR band. Handle is hidden so the legend contains curves only.
band = @(x,lo,hi,c) fill([x;flipud(x)], [lo;flipud(hi)], c, ...
    'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');

% -------------------------------------------------------------------------
% (a) Equilibrium composition error
% -------------------------------------------------------------------------
ax1 = nexttile; hold(ax1,'on');
band(T_grid,q25_pi(:,1),q75_pi(:,1),c1);
band(T_grid,q25_pi(:,2),q75_pi(:,2),c2);
h1 = plot(T_grid,med_pi(:,1),'-o','LineWidth',lw,'MarkerSize',ms,'Color',c1);
h2 = plot(T_grid,med_pi(:,2),'--s','LineWidth',lw,'MarkerSize',ms,'Color',c2);
xlabel('Observation length, $T$','Interpreter','latex');
ylabel('Median error in $\pi$','Interpreter','latex');
grid on; box on;
set(gca,'XScale','log','XTick',T_grid,'XTickLabel',string(T_grid));

% -------------------------------------------------------------------------
% (b) Mobility relative error
% -------------------------------------------------------------------------
ax2 = nexttile; hold(ax2,'on');
band(T_grid,q25_mu(:,1),q75_mu(:,1),c1);
band(T_grid,q25_mu(:,2),q75_mu(:,2),c2);
plot(T_grid,med_mu(:,1),'-o','LineWidth',lw,'MarkerSize',ms,'Color',c1);
plot(T_grid,med_mu(:,2),'--s','LineWidth',lw,'MarkerSize',ms,'Color',c2);
xlabel('Observation length, $T$','Interpreter','latex');
ylabel('Median relative error in $\mu$','Interpreter','latex');
grid on; box on;
set(gca,'XScale','log','XTick',T_grid,'XTickLabel',string(T_grid));

% -------------------------------------------------------------------------
% (c) Transition-matrix error
% -------------------------------------------------------------------------
ax3 = nexttile; hold(ax3,'on');
band(T_grid,q25_P(:,1),q75_P(:,1),c1);
band(T_grid,q25_P(:,2),q75_P(:,2),c2);
plot(T_grid,med_P(:,1),'-o','LineWidth',lw,'MarkerSize',ms,'Color',c1);
plot(T_grid,med_P(:,2),'--s','LineWidth',lw,'MarkerSize',ms,'Color',c2);
xlabel('Observation length, $T$','Interpreter','latex');
ylabel('Median error in $P$','Interpreter','latex');
grid on; box on;
set(gca,'XScale','log','XTick',T_grid,'XTickLabel',string(T_grid));

%% ---------------- SHARED LEGEND -----------------------------------------
lgd = legend(ax1,[h1 h2],graph_labels, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
try
    lgd.Layout.Tile = 'north';
catch
end

%% ---------------- TYPOGRAPHY --------------------------------------------
set([ax1 ax2 ax3], ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'TickLabelInterpreter','latex');

% Match the manuscript figure: major grid only.
for ax = [ax1 ax2 ax3]
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    xlim(ax,[min(T_grid) max(T_grid)]);
end

% Panel labels, positioned as in the T-sensitivity manuscript figure.
text(ax1,0.90,0.96,'(a)','Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',14);
text(ax2,0.90,0.96,'(b)','Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',14);
text(ax3,0.90,0.96,'(c)','Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',14);

% The shaded regions are the empirical 25th--75th percentile ranges across
% the R Monte Carlo replications. They represent dataset-to-dataset
% variability, not uncertainty in the median.

%% ---------------- PRINT NUMERICAL SUMMARY -------------------------------
fprintf('\nNonlocal-tree robustness: medians [IQR]\n');
fprintf('---------------------------------------------------------------\n');

for it = 1:numel(T_grid)
    fprintf('T = %d\n',T_grid(it));

    fprintf('  E_pi  Nearest-neighbour %.4f [%.4f, %.4f]   Nonlocal %.4f [%.4f, %.4f]\n', ...
        med_pi(it,1),q25_pi(it,1),q75_pi(it,1), ...
        med_pi(it,2),q25_pi(it,2),q75_pi(it,2));

    fprintf('  E_mu  Nearest-neighbour %.4f [%.4f, %.4f]   Nonlocal %.4f [%.4f, %.4f]\n', ...
        med_mu(it,1),q25_mu(it,1),q75_mu(it,1), ...
        med_mu(it,2),q25_mu(it,2),q75_mu(it,2));

    fprintf('  E_P   Nearest-neighbour %.4f [%.4f, %.4f]   Nonlocal %.4f [%.4f, %.4f]\n', ...
        med_P(it,1),q25_P(it,1),q75_P(it,1), ...
        med_P(it,2),q25_P(it,2),q75_P(it,2));
end

%% ---------------- SAVE --------------------------------------------------
savefig(fig,fig_file);

exportgraphics(fig,pdf_file,'ContentType','vector');
exportgraphics(fig,png_file,'Resolution',300);

fprintf('\nSaved:\n');
fprintf('  %s\n',fig_file);
fprintf('  %s\n',pdf_file);
fprintf('  %s\n',png_file);


%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function tf = is_connected_tree(E,K)
A=false(K,K);
for r=1:size(E,1)
    A(E(r,1),E(r,2))=true; A(E(r,2),E(r,1))=true;
end
seen=false(1,K); stack=1; seen(1)=true;
while ~isempty(stack)
    u=stack(end); stack(end)=[];
    nb=find(A(u,:));
    new=nb(~seen(nb)); seen(new)=true; stack=[stack new]; %#ok<AGROW>
end
tf=all(seen) && size(E,1)==K-1;
end

function v = normalize_positive(v)
v=max(v,realmin); v=v/sum(v);
end

function Q = build_generator_graph(pi_vec,mu,E)
pi_vec=pi_vec(:).'; mu=mu(:).';
K=numel(pi_vec); M=size(E,1);
if numel(mu)~=M, error('mu dimension must equal number of graph edges.'); end
Q=zeros(K,K);
for e=1:M
    u=E(e,1); v=E(e,2);
    a=mu(e)*sqrt(pi_vec(v)/pi_vec(u));
    b=mu(e)*sqrt(pi_vec(u)/pi_vec(v));
    Q(u,v)=Q(u,v)+a;
    Q(v,u)=Q(v,u)+b;
end
Q(1:K+1:end)=-sum(Q,2);
end

function P = transition_matrix_graph(pi_vec,mu,E,Delta)
P=expm(Delta*build_generator_graph(pi_vec,mu,E));
P(P<0 & P>-1e-13)=0; P=P./sum(P,2);
end

function theta = pack_theta_graph(pi_vec,mu)
K=numel(pi_vec);
theta=[log(pi_vec(1:K-1)./pi_vec(K)), log(mu(:).')].';
end

function [pi_vec,mu] = unpack_theta_graph(theta,K,M)
theta=theta(:);
if numel(theta)~=(K-1)+M, error('theta has wrong dimension.'); end
eta=[theta(1:K-1);0]; eta=eta-max(eta); w=exp(eta);
pi_vec=(w/sum(w)).';
mu=exp(theta(K:K+M-1)).';
end

function y = make_source_histogram_local(N,p)
a=N*p(:).'; y=floor(a); left=N-sum(y);
[~,ord]=sort(a-y,'descend');
if left>0, y(ord(1:left))=y(ord(1:left))+1; end
end

function Y = simulate_histogram_series_fast(y0,P,T)
y0=y0(:).'; K=numel(y0); Y=zeros(T+1,K); Y(1,:)=y0;
for t=1:T
    yn=zeros(1,K);
    for j=1:K
        if Y(t,j)>0, yn=yn+multinomial_sample_fast(Y(t,j),P(j,:)); end
    end
    Y(t+1,:)=yn;
end
end

function counts = multinomial_sample_fast(n,p)
p=max(p(:).',0); p=p/sum(p); K=numel(p); counts=zeros(1,K);
if n==0, return; end
if exist('binornd','file')==2
    remain_n=n; remain_p=1;
    for k=1:K-1
        if remain_n==0, break; end
        pk=min(max(p(k)/remain_p,0),1);
        counts(k)=binornd(remain_n,pk);
        remain_n=remain_n-counts(k); remain_p=remain_p-p(k);
        if remain_p<=eps, break; end
    end
    counts(K)=remain_n;
else
    edges=[0 cumsum(p)]; edges(end)=1;
    counts=histcounts(rand(n,1),edges);
end
end

function fit = fit_gaussian_graph(Y,E,Delta,start_theta,nstarts,maxit,maxfev)
K=size(Y,2); M=size(E,1);
opts=optimoptions('fminunc','Algorithm','quasi-newton', ...
    'SpecifyObjectiveGradient',true,'Display','off', ...
    'MaxIterations',maxit,'MaxFunctionEvaluations',maxfev, ...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-8,'FunctionTolerance',1e-8);
bestf=Inf; best=[];
for s=1:nstarts
    if s==1, th0=start_theta(:); else, th0=start_theta(:)+0.08*randn(size(start_theta(:))); end
    fun=@(th) gaussian_nll_score_graph(th,Y,E,Delta);
    [th,fval,ef,out]=fminunc(fun,th0,opts);
    if isfinite(fval) && fval<bestf
        bestf=fval; best.theta=th; best.exitflag=ef; best.output=out;
    end
end
if isempty(best), error('All Gaussian optimization starts failed.'); end
[pi_hat,mu_hat]=unpack_theta_graph(best.theta,K,M);
fit.theta=best.theta; fit.pi=pi_hat; fit.mu=mu_hat;
fit.P=transition_matrix_graph(pi_hat,mu_hat,E,Delta);
fit.loglik=-bestf; fit.exitflag=best.exitflag; fit.output=best.output;
end

function [nll,g] = gaussian_nll_score_graph(theta,Y,E,Delta)
K=size(Y,2); d=K-1; M=size(E,1); N=sum(Y(1,:));
[pi_vec,mu]=unpack_theta_graph(theta,K,M);
if any(~isfinite(pi_vec)) || any(~isfinite(mu)) || any(mu<=0) || max(mu)>1e6
    nll=1e100; g=zeros(size(theta)); return;
end
Q=build_generator_graph(pi_vec,mu,E); A=Delta*Q; P=expm(A);
if any(~isfinite(P(:))), nll=1e100; g=zeros(size(theta)); return; end
q=P(:,1:d); ll=0; GP=zeros(K,K);
for t=1:size(Y,1)-1
    rho=Y(t,:)/N; mfull=rho*P; m=mfull(1:d).';
    V=diag(m)-q.'*(rho(:).*q); V=(V+V.')/2;
    [R,pflag]=chol(V);
    if pflag~=0
        jitter=max(1e-12,1e-10*trace(V)/max(d,1)); V=V+jitter*eye(d); [R,pflag]=chol(V);
    end
    if pflag~=0, nll=1e100; g=zeros(size(theta)); return; end
    r=Y(t+1,1:d).'/N-m; z=R\(R'\r);
    ll=ll-0.5*(2*sum(log(diag(R)))+N*(r.'*z));
    Vinv=R\(R'\eye(d)); GV=0.5*(N*(z*z.')-Vinv); GV=(GV+GV.')/2;
    gm=N*z; dg=diag(GV);
    common=ones(K,1)*(gm+dg).'-2*q*GV;
    GP(:,1:d)=GP(:,1:d)+rho(:).*common;
end
Z=zeros(K,K); B=[A.',GP;Z,A.']; EB=expm(B);
GA=EB(1:K,K+1:2*K); GQ=Delta*GA;

% Chain Q -> log(pi) and edge-specific log(mu), generalized to arbitrary E.
g_s=zeros(K,1); g_logmu=zeros(M,1);
for e=1:M
    u=E(e,1); v=E(e,2);
    a=Q(u,v); b=Q(v,u);
    g_logmu(e)=a*(GQ(u,v)-GQ(u,u))+b*(GQ(v,u)-GQ(v,v));
    c=0.5*a*(GQ(u,u)-GQ(u,v))+0.5*b*(GQ(v,u)-GQ(v,v));
    g_s(u)=g_s(u)+c; g_s(v)=g_s(v)-c;
end
g_eta=g_s(1:K-1)-pi_vec(1:K-1).'*sum(g_s);
g_ll=[g_eta;g_logmu]; nll=-ll; g=-g_ll;
if ~isfinite(nll) || any(~isfinite(g)), nll=1e100; g=zeros(size(theta)); end
end

function check_gaussian_gradient_graph(theta,Y,E,Delta)
[f,g]=gaussian_nll_score_graph(theta,Y,E,Delta);
if ~isfinite(f), error('Gradient-check point has nonfinite objective.'); end
rng(9817,'twister'); m=numel(theta);
idx=unique(round(linspace(1,m,min(m,12)))); idx=unique([idx,randi(m,1,min(5,m))]);
rel=zeros(size(idx));
for ii=1:numel(idx)
    j=idx(ii); h=1e-6*(1+abs(theta(j))); ep=zeros(m,1); ep(j)=h;
    fp=gaussian_nll_score_graph(theta+ep,Y,E,Delta);
    fm=gaussian_nll_score_graph(theta-ep,Y,E,Delta);
    gn=(fp-fm)/(2*h);
    rel(ii)=abs(gn-g(j))/max([1,abs(gn),abs(g(j))]);
end
fprintf('  max relative coordinate-gradient discrepancy = %.3e\n',max(rel));
if max(rel)>5e-4, error('Analytic Gaussian gradient failed finite-difference check.'); end
end

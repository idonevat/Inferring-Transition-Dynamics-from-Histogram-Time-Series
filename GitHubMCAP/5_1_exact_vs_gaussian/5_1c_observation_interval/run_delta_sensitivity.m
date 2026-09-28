%% Experiment 1 diagnostic: sensitivity to observation interval Delta
%
% PURPOSE
% -------
% Test the interpretation that very large mobility estimates are primarily a
% temporal-resolution phenomenon.  The structural mechanism is held fixed and
% only the observation interval Delta is varied.
%
% Baseline model (same as the paper):
%   K  = 4
%   pi = (0.12, 0.28, 0.38, 0.22)
%   mu = (0.35, 0.25, 0.30)
%
% Diagnostic design:
%   N      = 20
%   T      = 20 observation intervals
%   Delta  = [0.25 0.5 1 2]
%   R      = 100 Monte Carlo replications
%
% Two estimators are fitted to every data set:
%   1. exact finite-population aggregate likelihood
%   2. Gaussian aggregate approximation
%
% Primary outputs:
%   - median relative Euclidean error in mu
%   - Pr(max_j mu_hat_j > 10)
%
% Secondary output:
%   - median relative Frobenius error in P = exp(Delta Q)
%
% IMPORTANT
% ---------
% This script is intentionally self-contained.  For speed, use the Parallel
% Computing Toolbox if available; set useParallel = true below.
%
% For a first check, set R = 20.  Once the code is verified, use R = 100 or
% R = 200 for the final diagnostic.
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% Reproducibility
masterSeed = 93457;
rng(masterSeed,'twister');

%% Design
K       = 4;
N       = 20;
T       = 20;
Delta   = [0.25 0.5 1 2];
R       = 100;

piTrue  = [0.12 0.28 0.38 0.22];
muTrue  = [0.35 0.25 0.30];

% Source-composition regime:
% 'equilibrium' is the cleanest diagnostic of temporal resolution.
sourceRegime = 'equilibrium';

% Definition of a fast-tail / effectively runaway estimate
muTailThreshold = 10;

% Matched mobility starts for both likelihoods
muStarts = [0.10 0.25 0.50];

% Optimizer controls
maxIter = 500;
maxFun  = 3000;
tolFun  = 1e-7;
tolX    = 1e-7;

% Parallel option
useParallel = false;

%% Derived truth
Qtrue = build_generator(piTrue,muTrue);

nD = numel(Delta);

% Results
muHatExact = nan(nD,R,K-1);
muHatGauss = nan(nD,R,K-1);

piHatExact = nan(nD,R,K);
piHatGauss = nan(nD,R,K);

PerrExact  = nan(nD,R);
PerrGauss  = nan(nD,R);

exitExact  = nan(nD,R);
exitGauss  = nan(nD,R);

timeExact  = nan(nD,R);
timeGauss  = nan(nD,R);

%% Main Monte Carlo loop
fprintf('\nDelta-sensitivity diagnostic\n');
fprintf('K=%d, N=%d, T=%d, R=%d\n',K,N,T,R);
fprintf('Delta = '); fprintf('%g ',Delta); fprintf('\n\n');

for d = 1:nD

    dt = Delta(d);
    Ptrue = expm(dt*Qtrue);

    fprintf('Delta = %.3g\n',dt);

    % Deterministic seeds make parfor/non-parfor runs reproducible.
    seeds = masterSeed + 100000*d + (1:R);

    % Temporary arrays for this Delta
    mE = nan(R,K-1);
    mG = nan(R,K-1);
    pE = nan(R,K);
    pG = nan(R,K);
    peE = nan(R,1);
    peG = nan(R,1);
    exE = nan(R,1);
    exG = nan(R,1);
    tmE = nan(R,1);
    tmG = nan(R,1);

    if useParallel

        parfor r = 1:R
            [mE(r,:),mG(r,:),pE(r,:),pG(r,:), ...
             peE(r),peG(r),exE(r),exG(r),tmE(r),tmG(r)] = ...
                one_replication(seeds(r),N,T,dt,piTrue,muTrue,Ptrue, ...
                                sourceRegime,muStarts, ...
                                maxIter,maxFun,tolFun,tolX);
        end

    else

        for r = 1:R
            [mE(r,:),mG(r,:),pE(r,:),pG(r,:), ...
             peE(r),peG(r),exE(r),exG(r),tmE(r),tmG(r)] = ...
                one_replication(seeds(r),N,T,dt,piTrue,muTrue,Ptrue, ...
                                sourceRegime,muStarts, ...
                                maxIter,maxFun,tolFun,tolX);

            if mod(r,10)==0 || r==R
                fprintf('  %4d / %4d complete\n',r,R);
            end
        end
    end

    muHatExact(d,:,:) = mE;
    muHatGauss(d,:,:) = mG;

    piHatExact(d,:,:) = pE;
    piHatGauss(d,:,:) = pG;

    PerrExact(d,:) = peE;
    PerrGauss(d,:) = peG;

    exitExact(d,:) = exE;
    exitGauss(d,:) = exG;

    timeExact(d,:) = tmE;
    timeGauss(d,:) = tmG;

    % Save after every Delta so a long run can be resumed/analyzed.
    save('experiment1_delta_sensitivity_checkpoint.mat', ...
         'K','N','T','Delta','R','piTrue','muTrue','Qtrue', ...
         'sourceRegime','muTailThreshold','muStarts', ...
         'muHatExact','muHatGauss','piHatExact','piHatGauss', ...
         'PerrExact','PerrGauss','exitExact','exitGauss', ...
         'timeExact','timeGauss','masterSeed','-v7.3');
end

%% Summaries
nmu = norm(muTrue);

muErrExact = nan(nD,R);
muErrGauss = nan(nD,R);

tailExact = false(nD,R);
tailGauss = false(nD,R);

for d = 1:nD
    AE = squeeze(muHatExact(d,:,:));  % R x (K-1)
    AG = squeeze(muHatGauss(d,:,:));

    muErrExact(d,:) = (vecnorm(AE-muTrue,2,2)/nmu).';
    muErrGauss(d,:) = (vecnorm(AG-muTrue,2,2)/nmu).';

    tailExact(d,:) = max(AE,[],2).' > muTailThreshold;
    tailGauss(d,:) = max(AG,[],2).' > muTailThreshold;
end

medMuExact = median(muErrExact,2,'omitnan').';
medMuGauss = median(muErrGauss,2,'omitnan').';

pctTailExact = 100*mean(tailExact,2,'omitnan').';
pctTailGauss = 100*mean(tailGauss,2,'omitnan').';

medPExact = median(PerrExact,2,'omitnan').';
medPGauss = median(PerrGauss,2,'omitnan').';

medTimeExact = median(timeExact,2,'omitnan').';
medTimeGauss = median(timeGauss,2,'omitnan').';

%% Print results
fprintf('\n============================================================\n');
fprintf('FINAL SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nDelta                  ');
fprintf('%10.3g',Delta);
fprintf('\n');

fprintf('Median rel mu err: exact');
fprintf('%10.4f',medMuExact);
fprintf('\n');

fprintf('Median rel mu err: gauss');
fprintf('%10.4f',medMuGauss);
fprintf('\n');

fprintf('Tail %% exact            ');
fprintf('%10.2f',pctTailExact);
fprintf('\n');

fprintf('Tail %% gaussian         ');
fprintf('%10.2f',pctTailGauss);
fprintf('\n');

fprintf('Median rel P err: exact ');
fprintf('%10.4f',medPExact);
fprintf('\n');

fprintf('Median rel P err: gauss ');
fprintf('%10.4f',medPGauss);
fprintf('\n');

fprintf('Median time exact (s)   ');
fprintf('%10.3f',medTimeExact);
fprintf('\n');

fprintf('Median time gauss (s)   ');
fprintf('%10.3f',medTimeGauss);
fprintf('\n');

%% Diagnostic figure
fig = figure('Color','w','Units','centimeters','Position',[2 2 17 7.5]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

lw = 1.5;
ms = 6;

% (a) mobility error
ax1 = nexttile(tl,1);
hold(ax1,'on');
plot(ax1,Delta,medMuExact,'-o','LineWidth',lw,'MarkerSize',ms, ...
    'DisplayName','Exact likelihood');
plot(ax1,Delta,medMuGauss,'-s','LineWidth',lw,'MarkerSize',ms, ...
    'DisplayName','Gaussian approximation');
grid(ax1,'on'); box(ax1,'on');
xlabel(ax1,'Observation interval, $\Delta$','Interpreter','latex');
ylabel(ax1,'Median relative error in $\mu$','Interpreter','latex');
title(ax1,'(a) Mobility estimation','FontWeight','normal');

% (b) fast-tail frequency
ax2 = nexttile(tl,2);
hold(ax2,'on');
plot(ax2,Delta,pctTailExact,'-o','LineWidth',lw,'MarkerSize',ms, ...
    'DisplayName','Exact likelihood');
plot(ax2,Delta,pctTailGauss,'-s','LineWidth',lw,'MarkerSize',ms, ...
    'DisplayName','Gaussian approximation');
grid(ax2,'on'); box(ax2,'on');
xlabel(ax2,'Observation interval, $\Delta$','Interpreter','latex');
ylabel(ax2,sprintf('Fits with $\\max_j\\hat{\\mu}_j>%g$ (\\%%)',muTailThreshold), ...
    'Interpreter','latex');
title(ax2,'(b) Fast-tail estimates','FontWeight','normal');

lgd = legend(ax1,'Location','northoutside','Orientation','horizontal');
lgd.Layout.Tile = 'north';

set([ax1 ax2], ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'LineWidth',0.75, ...
    'TickDir','out');

exportgraphics(fig,'experiment1_delta_sensitivity.pdf','ContentType','vector');
exportgraphics(fig,'experiment1_delta_sensitivity.png','Resolution',400);

%% Save final result
save('experiment1_delta_sensitivity_results.mat', ...
     'K','N','T','Delta','R','piTrue','muTrue','Qtrue', ...
     'sourceRegime','muTailThreshold','muStarts', ...
     'muHatExact','muHatGauss','piHatExact','piHatGauss', ...
     'PerrExact','PerrGauss','muErrExact','muErrGauss', ...
     'tailExact','tailGauss','exitExact','exitGauss', ...
     'timeExact','timeGauss', ...
     'medMuExact','medMuGauss','pctTailExact','pctTailGauss', ...
     'medPExact','medPGauss','medTimeExact','medTimeGauss', ...
     'masterSeed','-v7.3');

fprintf('\nSaved:\n');
fprintf('  experiment1_delta_sensitivity_results.mat\n');
fprintf('  experiment1_delta_sensitivity.pdf\n');
fprintf('  experiment1_delta_sensitivity.png\n');

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [muE,muG,piE,piG,pErrE,pErrG,exitE,exitG,tE,tG] = ...
    one_replication(seed,N,T,Delta,piTrue,muTrue,Ptrue, ...
                    sourceRegime,muStarts,maxIter,maxFun,tolFun,tolX)

    rng(seed,'twister');

    K = numel(piTrue);

    % ----- Initial aggregate histogram
    switch lower(sourceRegime)
        case 'equilibrium'
            rho0 = piTrue;

        case 'nonequilibrium'
            h = [0.12 -0.06 -0.04 -0.02];
            rho0 = piTrue + 0.5*h;
            rho0 = rho0/sum(rho0);

        otherwise
            error('Unknown sourceRegime: %s',sourceRegime);
    end

    y0 = integer_histogram(N,rho0);

    % ----- Simulate aggregate Markov trajectory exactly
    Y = simulate_histogram_chain(y0,T,Ptrue);

    % ----- pi starting value from empirical average composition
    piStart = mean(Y,1)/N;
    piStart = max(piStart,1e-4);
    piStart = piStart/sum(piStart);

    etaPi0 = log(piStart(1:K-1)/piStart(K));

    bestLE = -Inf;
    bestLG = -Inf;

    muE = nan(1,K-1);
    muG = nan(1,K-1);
    piE = nan(1,K);
    piG = nan(1,K);

    exitE = NaN;
    exitG = NaN;

    opts = optimoptions('fminunc', ...
        'Algorithm','quasi-newton', ...
        'Display','off', ...
        'MaxIterations',maxIter, ...
        'MaxFunctionEvaluations',maxFun, ...
        'OptimalityTolerance',tolFun, ...
        'StepTolerance',tolX, ...
        'FunctionTolerance',tolFun);

    % ----- Exact fit
    tic;
    for m0 = muStarts

        z0 = [etaPi0(:); log(m0)*ones(K-1,1)];

        obj = @(z) negloglik_exact(z,Y,Delta);

        try
            [zhat,fval,eflag] = fminunc(obj,z0,opts);

            if isfinite(fval) && -fval > bestLE
                bestLE = -fval;
                [piTmp,muTmp] = unpack_theta(zhat,K);
                piE = piTmp;
                muE = muTmp;
                exitE = eflag;
            end
        catch
            % Leave this start unused and continue to next matched start.
        end
    end
    tE = toc;

    % ----- Gaussian fit
    tic;
    for m0 = muStarts

        z0 = [etaPi0(:); log(m0)*ones(K-1,1)];

        obj = @(z) negloglik_gaussian(z,Y,Delta);

        try
            [zhat,fval,eflag] = fminunc(obj,z0,opts);

            if isfinite(fval) && -fval > bestLG
                bestLG = -fval;
                [piTmp,muTmp] = unpack_theta(zhat,K);
                piG = piTmp;
                muG = muTmp;
                exitG = eflag;
            end
        catch
        end
    end
    tG = toc;

    % ----- Relative P errors
    if all(isfinite(piE)) && all(isfinite(muE))
        QE = build_generator(piE,muE);
        PE = expm(Delta*QE);
        pErrE = norm(PE-Ptrue,'fro')/norm(Ptrue,'fro');
    else
        pErrE = NaN;
    end

    if all(isfinite(piG)) && all(isfinite(muG))
        QG = build_generator(piG,muG);
        PG = expm(Delta*QG);
        pErrG = norm(PG-Ptrue,'fro')/norm(Ptrue,'fro');
    else
        pErrG = NaN;
    end
end

% -------------------------------------------------------------------------
function Y = simulate_histogram_chain(y0,T,P)
% Exact simulation of the aggregate process without storing unit identities.
%
% Given current source counts y_j, each source group independently produces
% a multinomial destination vector with probabilities equal to row j of P.

    K = numel(y0);
    Y = zeros(T+1,K);
    Y(1,:) = y0;

    for t = 1:T
        y = Y(t,:);
        ynext = zeros(1,K);

        for j = 1:K
            if y(j)>0
                ynext = ynext + sample_multinomial(y(j),P(j,:));
            end
        end

        Y(t+1,:) = ynext;
    end
end

% -------------------------------------------------------------------------
function x = sample_multinomial(n,p)
% Toolbox-light multinomial sampler using sequential binomials.

    K = numel(p);
    x = zeros(1,K);

    remainingN = n;
    remainingP = 1;

    for k = 1:K-1
        if remainingN==0
            break;
        end

        pk = p(k)/remainingP;
        pk = min(max(pk,0),1);

        x(k) = binornd(remainingN,pk);

        remainingN = remainingN-x(k);
        remainingP = remainingP-p(k);

        if remainingP <= 0
            remainingP = eps;
        end
    end

    x(K) = remainingN;
end

% -------------------------------------------------------------------------
function y = integer_histogram(N,rho)
% Largest-remainder conversion of N*rho to a count histogram.

    a = N*rho;
    y = floor(a);
    remN = N-sum(y);

    frac = a-y;
    [~,idx] = sort(frac,'descend');

    if remN>0
        y(idx(1:remN)) = y(idx(1:remN))+1;
    end
end

% -------------------------------------------------------------------------
function nll = negloglik_exact(z,Y,Delta)

    K = size(Y,2);
    [pi,mu] = unpack_theta(z,K);

    Q = build_generator(pi,mu);
    P = expm(Delta*Q);

    % Reject numerical departures from a proper stochastic matrix.
    if any(~isfinite(P(:))) || any(P(:)<=0)
        nll = 1e100;
        return;
    end

    ll = 0;

    for t = 2:size(Y,1)
        y  = Y(t-1,:);
        yp = Y(t,:);

        pr = aggregate_transition_prob_dp(y,yp,P);

        if ~(isfinite(pr) && pr>0)
            nll = 1e100;
            return;
        end

        ll = ll + log(pr);
    end

    nll = -ll;
end

% -------------------------------------------------------------------------
function pr = aggregate_transition_prob_dp(y,target,P)
% Exact finite-population aggregate transition probability.
%
% The destination histogram is a sum of independent multinomial source-row
% contributions.  We convolve these row distributions successively while
% retaining only partial destination counts that do not exceed target.
%
% K=4 and N=20 make this diagnostic implementation practical.

    K = numel(y);

    % Encode a partial K-vector c in a scalar mixed-radix key.
    radix = target + 1;
    mult = ones(1,K);
    for k = 2:K
        mult(k) = mult(k-1)*radix(k-1);
    end

    zeroKey = 1;  % encoded zero vector +1 for MATLAB indexing
    maxKey  = prod(radix);

    dp = zeros(maxKey,1);
    dp(zeroKey) = 1;

    active = zeroKey;

    for j = 1:K

        nj = y(j);
        if nj==0
            continue;
        end

        rows = bounded_compositions(nj,target);
        nRows = size(rows,1);

        rowProb = zeros(nRows,1);
        for a = 1:nRows
            x = rows(a,:);
            rowProb(a) = multinomial_prob(x,P(j,:));
        end

        newdp = zeros(maxKey,1);
        newActiveFlag = false(maxKey,1);

        for aa = 1:numel(active)

            key = active(aa);
            base = decode_key(key-1,radix);

            pbase = dp(key);
            if pbase==0
                continue;
            end

            for b = 1:nRows

                c = base + rows(b,:);
                if any(c > target)
                    continue;
                end

                newKey = 1 + sum(c.*mult);
                newdp(newKey) = newdp(newKey) + pbase*rowProb(b);
                newActiveFlag(newKey) = true;
            end
        end

        dp = newdp;
        active = find(newActiveFlag);
    end

    targetKey = 1 + sum(target.*mult);
    pr = dp(targetKey);
end

% -------------------------------------------------------------------------
function C = bounded_compositions(n,upper)
% All K-part nonnegative integer compositions of n satisfying x<=upper.
% Recursive enumeration is adequate for K=4,N=20.

    K = numel(upper);
    buffer = zeros(10000,K);
    count = 0;
    x = zeros(1,K);

    recurse(1,n);

    C = buffer(1:count,:);

    function recurse(k,remaining)

        if k==K
            if remaining <= upper(k)
                x(k) = remaining;
                count = count+1;

                if count > size(buffer,1)
                    buffer = [buffer; zeros(size(buffer,1),K)]; %#ok<AGROW>
                end

                buffer(count,:) = x;
            end
            return;
        end

        maxHere = min(remaining,upper(k));

        for v = 0:maxHere
            x(k) = v;
            recurse(k+1,remaining-v);
        end
    end
end

% -------------------------------------------------------------------------
function p = multinomial_prob(x,prob)

    n = sum(x);

    logp = gammaln(n+1) - sum(gammaln(x+1));

    idx = x>0;
    if any(prob(idx)<=0)
        p = 0;
        return;
    end

    logp = logp + sum(x(idx).*log(prob(idx)));
    p = exp(logp);
end

% -------------------------------------------------------------------------
function c = decode_key(key0,radix)

    K = numel(radix);
    c = zeros(1,K);

    q = key0;

    for k = 1:K
        c(k) = mod(q,radix(k));
        q = floor(q/radix(k));
    end
end

% -------------------------------------------------------------------------
function nll = negloglik_gaussian(z,Y,Delta)

    K = size(Y,2);
    d = K-1;
    N = sum(Y(1,:));

    [pi,mu] = unpack_theta(z,K);
    Q = build_generator(pi,mu);
    P = expm(Delta*Q);

    if any(~isfinite(P(:))) || any(P(:)<=0)
        nll = 1e100;
        return;
    end

    ll = 0;

    for t = 2:size(Y,1)

        pprev = Y(t-1,:)/N;
        pnow  = Y(t,:)/N;

        mfull = pprev*P;
        m = mfull(1:d).';

        V = zeros(d,d);

        for j = 1:K
            q = P(j,1:d).';
            V = V + pprev(j)*(diag(q)-q*q.');
        end

        % Numerical safeguard
        V = (V+V.')/2;

        [L,pflag] = chol(V,'lower');

        if pflag~=0
            nll = 1e100;
            return;
        end

        r = pnow(1:d).' - m;

        % log det(V) via Cholesky
        logdetV = 2*sum(log(diag(L)));

        % r' V^{-1} r
        zres = L\r;
        quad = zres.'*zres;

        % Constants independent of theta omitted.
        ll = ll - 0.5*(logdetV + N*quad);
    end

    nll = -ll;
end

% -------------------------------------------------------------------------
function [pi,mu] = unpack_theta(z,K)
% Unconstrained coordinates:
%   pi_j / pi_K = exp(eta_j), j=1,...,K-1
%   mu_j = exp(lambda_j)

    eta = z(1:K-1);
    lam = z(K:end);

    % Stable softmax with last logit fixed to zero
    logits = [eta(:); 0];
    a = max(logits);
    w = exp(logits-a);
    pi = (w/sum(w)).';

    % Avoid floating-point overflow in pathological optimization proposals.
    lam = min(lam,log(1e14));
    mu = exp(lam).';
end

% -------------------------------------------------------------------------
function Q = build_generator(pi,mu)

    K = numel(pi);
    Q = zeros(K,K);

    for j = 1:K-1
        Q(j,j+1) = mu(j)*sqrt(pi(j+1)/pi(j));
        Q(j+1,j) = mu(j)*sqrt(pi(j)/pi(j+1));
    end

    for j = 1:K
        Q(j,j) = -sum(Q(j,:));
    end
end

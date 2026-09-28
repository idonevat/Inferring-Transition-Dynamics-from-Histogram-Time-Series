%EXPERIMENT1_PARAMETER_ROBUSTNESS_EXACT_VS_GAUSSIAN
% Robustness extension to Experiment 1: exact finite-population MLE versus
% Gaussian-approximate MLE across randomly generated transition mechanisms.
%
% REQUIRED helper functions from the existing Experiment 1 code folder:
%   transition_matrix.m, make_source_histogram.m,
%   simulate_histogram_series.m, pack_theta.m, fit_model.m,
%   prepare_exact_series_cache.m, and their dependencies.
%
% DESIGN
%   K = 4, Delta = 1, T = 20
%   Rpar independently drawn parameter mechanisms (pi,mu)
%   N in {5,20,40}
%   Balanced source regimes:
%       odd mechanism index  -> equilibrium source rho = pi
%       even mechanism index -> moderate nonequilibrium source
%   For every generated data set, fit BOTH exact and Gaussian MLEs using
%   identical starting values and identical random multistarts.
%
% The same randomly drawn mechanism is used across the N values. This makes
% comparisons across N less sensitive to changes in the mechanism mixture.
%
% The script checkpoints after every fitted data set and can be resumed.

clear; clc;
rng(20260908,'twister');
tic
%% Configuration
K = 4;
Delta = 1;
T = 20;
N_grid = [5 20 40];
Rpar = 100;             % number of independently drawn mechanisms
nstarts = 2;

% Parameter-generating distribution.
% pi is drawn from a symmetric Dirichlet(alpha_pi) distribution, with an
% interior rejection rule. alpha_pi is an integer so no Statistics Toolbox
% is needed by the helper below.
alpha_pi = 4;
pi_min = 0.05;
pi_max = 0.70;

% Mobilities are drawn independently log-uniformly over a moderate range.
% This spans slower/faster mechanisms around the paper's baseline values
% without intentionally concentrating on extreme ill-conditioning.
mu_min = 0.12;
mu_max = 0.60;

% For the moderate-nonequilibrium source, mix pi with an independent
% interior composition. lambda_source controls the displacement from pi.
lambda_source = 0.35;

checkpoint_file = 'experiment1_parameter_robustness_checkpoint.mat';
result_file = 'experiment1_parameter_robustness_results.mat';

%% Generate or restore the fixed mechanism bank
nN = numel(N_grid);

if exist(checkpoint_file,'file')
    S = load(checkpoint_file);
    required = {'raw','completed','pi_bank','mu_bank','rho_bank','regime_bank', ...
        'N_grid','T','Rpar','Delta','nstarts','alpha_pi','pi_min','pi_max', ...
        'mu_min','mu_max','lambda_source'};
    for jj = 1:numel(required)
        if ~isfield(S,required{jj})
            error('Checkpoint is missing field %s. Rename/delete it and restart.',required{jj});
        end
    end
    if ~isequal(S.N_grid,N_grid) || S.T~=T || S.Rpar~=Rpar || ...
            S.Delta~=Delta || S.nstarts~=nstarts || S.alpha_pi~=alpha_pi || ...
            S.pi_min~=pi_min || S.pi_max~=pi_max || S.mu_min~=mu_min || ...
            S.mu_max~=mu_max || S.lambda_source~=lambda_source
        error(['Existing robustness checkpoint uses a different configuration. ' ...
               'Rename/delete it before starting this design.']);
    end
    raw = S.raw;
    completed = S.completed;
    pi_bank = S.pi_bank;
    mu_bank = S.mu_bank;
    rho_bank = S.rho_bank;
    regime_bank = S.regime_bank;
    fprintf('Resuming checkpoint: %d/%d data sets already fitted.\n', ...
        nnz(completed),numel(completed));
else
    pi_bank = nan(Rpar,K);
    mu_bank = nan(Rpar,K-1);
    rho_bank = nan(Rpar,K);
    regime_bank = strings(Rpar,1);

    for r = 1:Rpar
        pi_r = draw_interior_dirichlet(K,alpha_pi,pi_min,pi_max);
        logmu = log(mu_min) + (log(mu_max)-log(mu_min))*rand(1,K-1);
        mu_r = exp(logmu);

        if mod(r,2)==1
            rho_r = pi_r;
            regime_bank(r) = "equilibrium";
        else
            rho_alt = draw_interior_dirichlet(K,alpha_pi,pi_min,pi_max);
            rho_r = (1-lambda_source)*pi_r + lambda_source*rho_alt;
            rho_r = rho_r/sum(rho_r);
            regime_bank(r) = "non-equilibrium";
        end

        pi_bank(r,:) = pi_r;
        mu_bank(r,:) = mu_r;
        rho_bank(r,:) = rho_r;
    end

    raw = initialize_raw(Rpar,nN,K);
    completed = false(Rpar,nN);
end

%% Monte Carlo over mechanisms and population sizes
for r = 1:Rpar
    pi_true = pi_bank(r,:);
    mu_true = mu_bank(r,:);
    rho = rho_bank(r,:);
    P_true = transition_matrix(pi_true,mu_true,Delta);

    for iN = 1:nN
        if completed(r,iN)
            continue;
        end

        N = N_grid(iN);
        y0 = make_source_histogram(N,rho);
        Y = simulate_histogram_series(y0,P_true,T);

        % Common mildly perturbed start around truth. This experiment studies
        % estimator/likelihood behavior rather than initialization robustness.
        pi_start = pi_true.*exp(0.05*randn(1,K));
        pi_start = pi_start/sum(pi_start);
        mu_start = mu_true.*exp(0.05*randn(1,K-1));
        start = pack_theta(pi_start,mu_start);

        % Build exact combinatorial cache once for the observed series.
        tic;
        exact_cache = prepare_exact_series_cache(Y);
        raw.cache_time(r,iN) = toc;

        % Match the random multistarts exactly across the two estimators.
        multistart_state = rng;
        tic;
        fE = fit_model(Y,Delta,'exact_cached',start,nstarts,exact_cache);
        raw.opt_time_exact(r,iN) = toc;

        rng(multistart_state);
        tic;
        fG = fit_model(Y,Delta,'gaussian',start,nstarts);
        raw.opt_time_gaussian(r,iN) = toc;

        raw.end_time_exact(r,iN) = raw.cache_time(r,iN) + raw.opt_time_exact(r,iN);

        % Relative errors to truth. These are directly interpretable across
        % mechanisms whose parameter magnitudes differ.
        raw.rel_pi_exact(r,iN) = norm(fE.pi-pi_true)/norm(pi_true);
        raw.rel_mu_exact(r,iN) = norm(fE.mu-mu_true)/norm(mu_true);
        raw.rel_P_exact(r,iN) = norm(fE.P-P_true,'fro')/norm(P_true,'fro');

        raw.rel_pi_gaussian(r,iN) = norm(fG.pi-pi_true)/norm(pi_true);
        raw.rel_mu_gaussian(r,iN) = norm(fG.mu-mu_true)/norm(mu_true);
        raw.rel_P_gaussian(r,iN) = norm(fG.P-P_true,'fro')/norm(P_true,'fro');

        % Direct GA-versus-exact discrepancy. Normalize by the exact fitted
        % quantity, with a small numerical floor only to guard division.
        raw.direct_pi(r,iN) = norm(fG.pi-fE.pi)/max(norm(fE.pi),1e-12);
        raw.direct_mu(r,iN) = norm(fG.mu-fE.mu)/max(norm(fE.mu),1e-12);
        raw.direct_P(r,iN) = norm(fG.P-fE.P,'fro')/max(norm(fE.P,'fro'),1e-12);

        % Signed difference in error to truth: positive means GA has larger
        % relative error than exact MLE on this data set.
        raw.penalty_pi(r,iN) = raw.rel_pi_gaussian(r,iN)-raw.rel_pi_exact(r,iN);
        raw.penalty_mu(r,iN) = raw.rel_mu_gaussian(r,iN)-raw.rel_mu_exact(r,iN);
        raw.penalty_P(r,iN) = raw.rel_P_gaussian(r,iN)-raw.rel_P_exact(r,iN);

        raw.exit_exact(r,iN) = fE.exitflag;
        raw.exit_gaussian(r,iN) = fG.exitflag;
        raw.max_mu_exact(r,iN) = max(fE.mu);
        raw.max_mu_gaussian(r,iN) = max(fG.mu);

        completed(r,iN) = true;

        save(checkpoint_file,'raw','completed','pi_bank','mu_bank','rho_bank', ...
            'regime_bank','N_grid','T','Rpar','Delta','nstarts','alpha_pi', ...
            'pi_min','pi_max','mu_min','mu_max','lambda_source','-v7.3');

        fprintf(['mechanism=%3d/%3d %-15s N=%2d | exact=%7.2fs gauss=%6.2fs ' ...
                 '| rel-mu=(%.3f, %.3f) direct-mu=%.3f exit=(%d,%d)\n'], ...
            r,Rpar,char(regime_bank(r)),N,raw.end_time_exact(r,iN), ...
            raw.opt_time_gaussian(r,iN),raw.rel_mu_exact(r,iN), ...
            raw.rel_mu_gaussian(r,iN),raw.direct_mu(r,iN), ...
            fE.exitflag,fG.exitflag);
    end
end

%% Summary tables
summary_all = summarize_subset(raw,true(Rpar,1),N_grid,'all');
summary_eq = summarize_subset(raw,regime_bank=="equilibrium",N_grid,'equilibrium');
summary_neq = summarize_subset(raw,regime_bank=="non-equilibrium",N_grid,'non-equilibrium');

save(result_file,'raw','summary_all','summary_eq','summary_neq','pi_bank','mu_bank', ...
    'rho_bank','regime_bank','N_grid','T','Rpar','Delta','nstarts','alpha_pi', ...
    'pi_min','pi_max','mu_min','mu_max','lambda_source','-v7.3');

fprintf('\n===============================================================\n');
fprintf('PARAMETER-ROBUSTNESS STUDY COMPLETE\n');
fprintf('===============================================================\n');
disp(summary_all);
fprintf('\nEquilibrium-source subset:\n');
disp(summary_eq);
fprintf('\nNon-equilibrium-source subset:\n');
disp(summary_neq);
fprintf('Results saved to %s\n',result_file);

%% Simple diagnostic figure
% Three panels: mobility error, P error, and direct P discrepancy.
figure('Color','w','Position',[100 100 1200 350]);

subplot(1,3,1); hold on;
medE = median(raw.rel_mu_exact,1,'omitnan');
medG = median(raw.rel_mu_gaussian,1,'omitnan');
plot(N_grid,medE,'-o','LineWidth',1.5,'DisplayName','Exact MLE');
plot(N_grid,medG,'--s','LineWidth',1.5,'DisplayName','Gaussian MLE');
xlabel('N'); ylabel('Median relative error in \mu'); grid on; legend('Location','best');
title('Transition intensity');

subplot(1,3,2); hold on;
medE = median(raw.rel_P_exact,1,'omitnan');
medG = median(raw.rel_P_gaussian,1,'omitnan');
plot(N_grid,medE,'-o','LineWidth',1.5,'DisplayName','Exact MLE');
plot(N_grid,medG,'--s','LineWidth',1.5,'DisplayName','Gaussian MLE');
xlabel('N'); ylabel('Median relative error in P'); grid on; legend('Location','best');
title('Finite-time transition matrix');

subplot(1,3,3);
boxplot(raw.direct_P,'Labels',string(N_grid));
xlabel('N'); ylabel('Relative distance: GA versus exact P'); grid on;
title('Approximation discrepancy');

sgtitle(sprintf('Robustness across %d randomly generated mechanisms, T=%d',Rpar,T));
saveas(gcf,'experiment1_parameter_robustness_diagnostic.png');
toc
%% ================================================================
% Local helpers
% ================================================================
function p = draw_interior_dirichlet(K,alpha,pmin,pmax)
% Draw Dirichlet(alpha,...,alpha) without requiring Statistics Toolbox.
% alpha must be a positive integer here. A Gamma(alpha,1) variate is the
% sum of alpha independent Exp(1) variables.
if alpha~=round(alpha) || alpha<1
    error('alpha must be a positive integer in this toolbox-free helper.');
end
for attempt = 1:100000
    g = zeros(1,K);
    for a = 1:alpha
        g = g - log(max(rand(1,K),realmin));
    end
    p = g/sum(g);
    if min(p)>=pmin && max(p)<=pmax
        return;
    end
end
error('Could not draw an interior Dirichlet composition; relax pmin/pmax.');
end

function raw = initialize_raw(Rpar,nN,K)
fields = {'rel_pi_exact','rel_mu_exact','rel_P_exact', ...
          'rel_pi_gaussian','rel_mu_gaussian','rel_P_gaussian', ...
          'direct_pi','direct_mu','direct_P', ...
          'penalty_pi','penalty_mu','penalty_P', ...
          'cache_time','opt_time_exact','end_time_exact','opt_time_gaussian', ...
          'exit_exact','exit_gaussian','max_mu_exact','max_mu_gaussian'};
for j=1:numel(fields)
    raw.(fields{j}) = nan(Rpar,nN);
end
raw.K = K;
end

function S = summarize_subset(raw,idx,N_grid,label)
nN = numel(N_grid);
S = table(N_grid(:),'VariableNames',{'N'});
S.subset = repmat(string(label),nN,1);
S.n = repmat(sum(idx),nN,1);

names = {'rel_pi_exact','rel_pi_gaussian','rel_mu_exact','rel_mu_gaussian', ...
         'rel_P_exact','rel_P_gaussian','direct_mu','direct_P', ...
         'penalty_mu','penalty_P','end_time_exact','opt_time_gaussian'};
for j=1:numel(names)
    A = raw.(names{j})(idx,:);
    S.(['median_' names{j}]) = median(A,1,'omitnan')';
end

S.q90_rel_mu_exact = percentile_cols(raw.rel_mu_exact(idx,:),90)';
S.q90_rel_mu_gaussian = percentile_cols(raw.rel_mu_gaussian(idx,:),90)';
S.q90_rel_P_exact = percentile_cols(raw.rel_P_exact(idx,:),90)';
S.q90_rel_P_gaussian = percentile_cols(raw.rel_P_gaussian(idx,:),90)';
S.q90_direct_mu = percentile_cols(raw.direct_mu(idx,:),90)';
S.q90_direct_P = percentile_cols(raw.direct_P(idx,:),90)';
S.exact_success_rate = mean(raw.exit_exact(idx,:)>0,1,'omitnan')';
S.gaussian_success_rate = mean(raw.exit_gaussian(idx,:)>0,1,'omitnan')';
S.frac_mu_exact_gt10 = mean(raw.max_mu_exact(idx,:)>10,1,'omitnan')';
S.frac_mu_gaussian_gt10 = mean(raw.max_mu_gaussian(idx,:)>10,1,'omitnan')';
S.speedup_exact_over_gaussian = S.median_end_time_exact ./ S.median_opt_time_gaussian;
end

function q = percentile_cols(A,p)
% Toolbox-free percentile by columns using linear interpolation.
q = nan(1,size(A,2));
for c=1:size(A,2)
    x=A(:,c); x=x(isfinite(x)); x=sort(x);
    n=numel(x);
    if n==0, continue; end
    if n==1, q(c)=x; continue; end
    pos=1+(n-1)*p/100;
    lo=floor(pos); hi=ceil(pos);
    if lo==hi
        q(c)=x(lo);
    else
        q(c)=x(lo)+(pos-lo)*(x(hi)-x(lo));
    end
end
end

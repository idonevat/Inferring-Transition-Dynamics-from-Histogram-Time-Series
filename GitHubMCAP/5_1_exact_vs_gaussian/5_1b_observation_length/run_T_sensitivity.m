%RUN_T_SENSITIVITY Effect of observation length T on exact and Gaussian MLE.
%
% This is an add-on to Experiment 1. It keeps population size fixed and
% varies the number of observed histogram transitions. Both exact cached-DP
% and Gaussian likelihoods estimate (pi,mu) using the same parameterization,
% starting point, and random multistarts.
%
% Design:
%   K = 4, Delta = 1, N = 20
%   T in {10,20,50,100,200}
%   R = 100 paired Monte Carlo trajectories per source regime
%   source regimes: equilibrium and moderate nonequilibrium
%
% Each replication simulates one trajectory of length max(T_grid), then fits
% every prefix. This gives a paired comparison across T and avoids changing
% the underlying realization when T is increased.
%
% Requires the functions in ../publication_support (or that folder on path).

clear; clc;
rng(2031,'twister');

%% Locate shared publication functions
this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
support_candidates = { ...
    fullfile(this_dir,'..','publication_support'), ...
    fullfile(this_dir,'..','code','publication_support'), ...
    fullfile(this_dir,'..','..','publication_support'), ...
    fullfile(this_dir,'..','..','code','publication_support')};
for ii=1:numel(support_candidates)
    if isfolder(support_candidates{ii})
        addpath(support_candidates{ii});
        break;
    end
end
assert(exist('fit_model','file')==2, ...
    'Could not find publication_support. Add that folder to the MATLAB path.');

%% Configuration
K = 4;
Delta = 1;
N = 20;
T_grid = [10 20 50 100 200];
Tmax = max(T_grid);
R = 100;
nstarts = 2;

% Baseline mechanism stated in the manuscript.
pi_true = [0.12 0.28 0.38 0.22];
mu_true = [0.35 0.25 0.30];
P_true = transition_matrix(pi_true,mu_true,Delta);

h = [0.12 -0.06 -0.04 -0.02];
delta_grid = [0 0.5];
scenario_name = {'equilibrium','non-equilibrium'};

checkpoint_file = 'experiment1_T_sensitivity_checkpoint.mat';
result_file = 'experiment1_T_sensitivity_results.mat';

nS = numel(delta_grid);
nT = numel(T_grid);

%% Allocate or resume
if isfile(checkpoint_file)
    S = load(checkpoint_file);
    required = {'raw','completed','T_grid','N','R','delta_grid','pi_true', ...
        'mu_true','Delta','nstarts'};
    for jj=1:numel(required)
        assert(isfield(S,required{jj}),'Checkpoint missing field %s.',required{jj});
    end
    assert(isequal(S.T_grid,T_grid) && S.N==N && S.R==R && ...
        isequal(S.delta_grid,delta_grid) && isequal(S.pi_true,pi_true) && ...
        isequal(S.mu_true,mu_true) && S.Delta==Delta && S.nstarts==nstarts, ...
        'Existing checkpoint uses a different configuration.');
    raw = S.raw;
    completed = S.completed;
    fprintf('Resuming: %d/%d fits complete.\n',nnz(completed),numel(completed));
else
    raw = initialize_raw(nS,nT,R,K);
    completed = false(nS,nT,R);
end

%% Monte Carlo loop
for s=1:nS
    rho = pi_true + delta_grid(s)*h;
    assert(all(rho>0),'Source composition is not interior.');
    rho = rho/sum(rho);

    for r=1:R
        % If every T for this scenario/replication is complete, skip simulation.
        if all(completed(s,:,r))
            continue;
        end

        rep_seed = 20310000 + 100000*s + r;
        rng(rep_seed,'twister');

        y0 = make_source_histogram(N,rho);
        Yfull = simulate_histogram_series(y0,P_true,Tmax);

        % One common base start for every T in this replication.
        pi_start = pi_true.*exp(0.05*randn(1,K));
        pi_start = pi_start/sum(pi_start);
        mu_start = mu_true.*exp(0.05*randn(1,K-1));
        start = pack_theta(pi_start,mu_start);

        % Save RNG state so every T and both likelihoods use identical
        % multistart perturbations within this replication.
        multistart_state = rng;

        for iT=1:nT
            if completed(s,iT,r)
                continue;
            end

            T = T_grid(iT);
            Y = Yfull(1:T+1,:);

            tic;
            exact_cache = prepare_exact_series_cache(Y);
            raw.cache_time(s,iT,r) = toc;

            rng(multistart_state);
            tic;
            fE = fit_model(Y,Delta,'exact_cached',start,nstarts,exact_cache);
            raw.opt_time_exact(s,iT,r) = toc;

            rng(multistart_state);
            tic;
            fG = fit_model(Y,Delta,'gaussian',start,nstarts);
            raw.opt_time_gaussian(s,iT,r) = toc;

            raw.end_time_exact(s,iT,r) = raw.cache_time(s,iT,r) + raw.opt_time_exact(s,iT,r);

            % Store estimates.
            raw.pi_hat_exact(s,iT,r,:) = reshape(fE.pi,1,1,1,K);
            raw.pi_hat_gaussian(s,iT,r,:) = reshape(fG.pi,1,1,1,K);
            raw.mu_hat_exact(s,iT,r,:) = reshape(fE.mu,1,1,1,K-1);
            raw.mu_hat_gaussian(s,iT,r,:) = reshape(fG.mu,1,1,1,K-1);

            % Errors use the definitions stated in the manuscript.
            raw.err_pi_exact(s,iT,r) = norm(fE.pi-pi_true);
            raw.err_pi_gaussian(s,iT,r) = norm(fG.pi-pi_true);
            raw.err_P_exact(s,iT,r) = norm(fE.P-P_true,'fro');
            raw.err_P_gaussian(s,iT,r) = norm(fG.P-P_true,'fro');
            raw.err_mu_exact(s,iT,r) = norm(fE.mu-mu_true)/norm(mu_true);
            raw.err_mu_gaussian(s,iT,r) = norm(fG.mu-mu_true)/norm(mu_true);

            % Direct GA-versus-exact discrepancy.
            raw.diff_pi(s,iT,r) = norm(fG.pi-fE.pi);
            raw.diff_P(s,iT,r) = norm(fG.P-fE.P,'fro');
            raw.diff_mu(s,iT,r) = norm(fG.mu-fE.mu)/norm(mu_true);

            raw.max_mu_exact(s,iT,r) = max(fE.mu);
            raw.max_mu_gaussian(s,iT,r) = max(fG.mu);
            raw.exit_exact(s,iT,r) = fE.exitflag;
            raw.exit_gaussian(s,iT,r) = fG.exitflag;

            completed(s,iT,r) = true;
            save(checkpoint_file,'raw','completed','T_grid','N','R','delta_grid', ...
                'scenario_name','pi_true','mu_true','P_true','Delta','nstarts','h','-v7.3');

            fprintf(['scenario=%-15s T=%3d rep=%3d/%3d | cache=%7.2fs ' ...
                'exact=%7.2fs gauss=%7.2fs | exit=(%d,%d)\n'], ...
                scenario_name{s},T,r,R,raw.cache_time(s,iT,r), ...
                raw.opt_time_exact(s,iT,r),raw.opt_time_gaussian(s,iT,r), ...
                fE.exitflag,fG.exitflag);
        end
    end
end

%% Summaries
results = struct();
results.scenario_name = scenario_name;
results.T_grid = T_grid;
results.N = N;

metric_names = {'err_pi','err_P','err_mu','diff_pi','diff_P','diff_mu'};
for s=1:nS
    for iT=1:nT
        z.scenario = scenario_name{s}; %#ok<AGROW>
        z.T = T_grid(iT);
        for mm=1:numel(metric_names)
            base = metric_names{mm};
            if startsWith(base,'err_')
                xE = squeeze(raw.([base '_exact'])(s,iT,:));
                xG = squeeze(raw.([base '_gaussian'])(s,iT,:));
                z.(['median_' base '_exact']) = median(xE,'omitnan');
                z.(['median_' base '_gaussian']) = median(xG,'omitnan');
                z.(['rmse_' base '_exact']) = sqrt(mean(xE.^2,'omitnan'));
                z.(['rmse_' base '_gaussian']) = sqrt(mean(xG.^2,'omitnan'));
            else
                x = squeeze(raw.(base)(s,iT,:));
                z.(['median_' base]) = median(x,'omitnan');
            end
        end
        z.frac_maxmu_gt10_exact = mean(squeeze(raw.max_mu_exact(s,iT,:))>10,'omitnan');
        z.frac_maxmu_gt10_gaussian = mean(squeeze(raw.max_mu_gaussian(s,iT,:))>10,'omitnan');
        z.exact_success_rate = mean(squeeze(raw.exit_exact(s,iT,:))>0,'omitnan');
        z.gaussian_success_rate = mean(squeeze(raw.exit_gaussian(s,iT,:))>0,'omitnan');
        z.median_time_exact = median(squeeze(raw.end_time_exact(s,iT,:)),'omitnan');
        z.median_time_gaussian = median(squeeze(raw.opt_time_gaussian(s,iT,:)),'omitnan');
        results.summary(s,iT) = z; %#ok<SAGROW>
        clear z
    end
end

save(result_file,'results','raw','pi_true','mu_true','P_true','T_grid','N', ...
    'delta_grid','scenario_name','R','Delta','nstarts','h','-v7.3');

fprintf('\nComplete. Saved %s\n',result_file);
fprintf('Run plot_T_sensitivity after inspecting the diagnostics.\n');

%% Local helper
function raw = initialize_raw(nS,nT,R,K)
raw.err_pi_exact = nan(nS,nT,R);
raw.err_pi_gaussian = nan(nS,nT,R);
raw.err_P_exact = nan(nS,nT,R);
raw.err_P_gaussian = nan(nS,nT,R);
raw.err_mu_exact = nan(nS,nT,R);
raw.err_mu_gaussian = nan(nS,nT,R);
raw.diff_pi = nan(nS,nT,R);
raw.diff_P = nan(nS,nT,R);
raw.diff_mu = nan(nS,nT,R);
raw.pi_hat_exact = nan(nS,nT,R,K);
raw.pi_hat_gaussian = nan(nS,nT,R,K);
raw.mu_hat_exact = nan(nS,nT,R,K-1);
raw.mu_hat_gaussian = nan(nS,nT,R,K-1);
raw.max_mu_exact = nan(nS,nT,R);
raw.max_mu_gaussian = nan(nS,nT,R);
raw.cache_time = nan(nS,nT,R);
raw.opt_time_exact = nan(nS,nT,R);
raw.end_time_exact = nan(nS,nT,R);
raw.opt_time_gaussian = nan(nS,nT,R);
raw.exit_exact = nan(nS,nT,R);
raw.exit_gaussian = nan(nS,nT,R);
end

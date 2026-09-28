%EXPERIMENT1_EXACT_VS_GAUSSIAN Exact cached-DP versus Gaussian MLE.
% Publication-design experiment with checkpoint/resume support.
%
% Main design:
%   Frozen reference mechanism 93 from the heterogeneous-mechanism bank.
%   It was selected using true model characteristics only, before examining
%   its estimation performance. N in {3,5,7,10,20,30,40}, T=20, two source regimes,
%   exact cached likelihood versus Gaussian approximation.
%
% The script saves experiment1_checkpoint.mat after EVERY replication.
% If an R=30 checkpoint exists, it is backed up and extended to R=100.
% If interrupted, rerunning resumes from the first unfinished replication.

clear; clc;
rng(2027,'twister');

%% Configuration
K = 4;
Delta = 1;
pi_true = [0.26660926 0.13909125 0.35642353 0.23787596];
mu_true = [0.25738153 0.46589071 0.18508825];
P_true = transition_matrix(pi_true,mu_true,Delta);

N_grid = [3 5 7 10 20 30 40];
T = 20;
R = 200;                % Pilot size; increase only after inspecting tails.
nstarts = 2;

% Two source regimes. The non-equilibrium perturbation is deliberately
% moderate: rho = pi + 0.5*h, with h*1 = 0.
h = [0.12 -0.06 -0.04 -0.02];
delta_grid = [0 0.5];
scenario_name = {'equilibrium','non-equilibrium'};

checkpoint_file = 'experiment1_reference_m93_R200_checkpoint.mat';
result_file = 'experiment1_reference_m93_R200_results.mat';

%% Allocate or resume
nS = numel(delta_grid);
nN = numel(N_grid);

if exist(checkpoint_file,'file')
    S = load(checkpoint_file);
    required = {'raw','completed','N_grid','T','R','delta_grid','pi_true','mu_true','Delta','nstarts'};
    for jj=1:numel(required)
        if ~isfield(S,required{jj})
            error('Checkpoint is missing field %s. Rename/delete it and restart.',required{jj});
        end
    end

    % Configuration must agree except that an older checkpoint may have
    % fewer Monte Carlo replications.
    if ~isequal(S.N_grid,N_grid) || S.T~=T || ...
            ~isequal(S.delta_grid,delta_grid) || ~isequal(S.pi_true,pi_true) || ...
            ~isequal(S.mu_true,mu_true) || S.Delta~=Delta || S.nstarts~=nstarts
        error(['Existing mechanism-93 checkpoint uses a different configuration. ' ...
               'Rename/delete it before starting this design.']);
    end

    if S.R > R
        error('Checkpoint has R=%d, which exceeds the requested R=%d.',S.R,R);
    end

    oldR = S.R;

    % Preserve a copy of the original smaller checkpoint before extending it.
    if oldR < R
        backup_file = sprintf('experiment1_reference_m93_checkpoint_R%d.mat',oldR);
        if ~isfile(backup_file)
            copyfile(checkpoint_file,backup_file);
            fprintf('Backed up original checkpoint as %s\n',backup_file);
        end

        raw_old = S.raw;
        completed_old = S.completed;

        raw = initialize_raw(nS,nN,R,K);

        numeric_fields = {'se_pi_exact','se_mu_exact','se_P_exact', ...
            'se_pi_gaussian','se_mu_gaussian','se_P_gaussian', ...
            'd2_pi','d2_mu','d2_P','cache_time','opt_time_exact', ...
            'end_time_exact','opt_time_gaussian','exit_exact','exit_gaussian', ...
            'mu_hat_exact','mu_hat_gaussian','max_mu_exact','max_mu_gaussian'};

        for jj=1:numel(numeric_fields)
            f = numeric_fields{jj};
            if isfield(raw_old,f)
                if ndims(raw.(f))==4
                    raw.(f)(:,:,1:oldR,:) = raw_old.(f);
                else
                    raw.(f)(:,:,1:oldR) = raw_old.(f);
                end
            end
        end

        raw.P_hat_exact = raw_old.P_hat_exact;
        raw.P_hat_gaussian = raw_old.P_hat_gaussian;
        raw.Y_example = raw_old.Y_example;

        completed = false(nS,nN,R);
        completed(:,:,1:oldR) = completed_old;

        fprintf('Extended checkpoint from R=%d to R=%d. Existing runs are preserved.\n',oldR,R);

        % Save the expanded checkpoint immediately.
        save(checkpoint_file,'raw','completed','N_grid','T','R','delta_grid', ...
            'scenario_name','pi_true','mu_true','P_true','Delta','nstarts','h','-v7.3');
    else
        raw = S.raw;
        completed = S.completed;
        oldR = R;
    end

    fprintf('Resuming checkpoint: %d/%d replications already complete.\n', ...
        nnz(completed),numel(completed));
else
    oldR = 0;
    raw = initialize_raw(nS,nN,R,K);
    completed = false(nS,nN,R);
end

%% Monte Carlo loop
for s = 1:nS
    delta = delta_grid(s);
    rho = pi_true + delta*h;
    if any(rho<=0)
        error('Scenario %d produces a non-positive source composition.',s);
    end
    rho = rho/sum(rho);

    for iN = 1:nN
        N = N_grid(iN);

        for r = 1:R
            if completed(s,iN,r)
                continue;
            end

            % For new replications, use a deterministic replication-specific
            % seed. This avoids reusing the original rng(2026) stream after
            % skipping the preserved R=30 runs.
            rep_seed = 20260000 + 100000*s + 1000*iN + r;
            rng(rep_seed,'twister');

            % Deterministic source histogram nearest N*rho, then exact aggregate
            % simulation under the same unit-level transition matrix.
            y0 = make_source_histogram(N,rho);
            Y = simulate_histogram_series(y0,P_true,T);

            % Mildly perturbed common first start.
            pi_start = pi_true.*exp(0.05*randn(1,K));
            pi_start = pi_start/sum(pi_start);
            mu_start = mu_true.*exp(0.05*randn(1,K-1));
            start = pack_theta(pi_start,mu_start);

            % Build exact combinatorial cache once for this observed series.
            tic;
            exact_cache = prepare_exact_series_cache(Y);
            raw.cache_time(s,iN,r) = toc;

            % Force exact and Gaussian fits to use the same random multistarts.
            multistart_state = rng;
            tic;
            fE = fit_model(Y,Delta,'exact_cached',start,nstarts,exact_cache);
            raw.opt_time_exact(s,iN,r) = toc;

            rng(multistart_state);
            tic;
            fG = fit_model(Y,Delta,'gaussian',start,nstarts);
            raw.opt_time_gaussian(s,iN,r) = toc;

            raw.end_time_exact(s,iN,r) = raw.cache_time(s,iN,r) + raw.opt_time_exact(s,iN,r);

            % Squared estimation errors, normalized by parameter dimension.
            raw.se_pi_exact(s,iN,r) = norm(fE.pi-pi_true)^2/K;
            raw.se_mu_exact(s,iN,r) = norm(fE.mu-mu_true)^2/(K-1);
            raw.se_P_exact(s,iN,r) = norm(fE.P-P_true,'fro')^2/K^2;

            raw.se_pi_gaussian(s,iN,r) = norm(fG.pi-pi_true)^2/K;
            raw.se_mu_gaussian(s,iN,r) = norm(fG.mu-mu_true)^2/(K-1);
            raw.se_P_gaussian(s,iN,r) = norm(fG.P-P_true,'fro')^2/K^2;

            % Direct exact-versus-Gaussian estimator disagreement.
            raw.d2_pi(s,iN,r) = norm(fE.pi-fG.pi)^2/K;
            raw.d2_mu(s,iN,r) = norm(fE.mu-fG.mu)^2/(K-1);
            raw.d2_P(s,iN,r) = norm(fE.P-fG.P,'fro')^2/K^2;

            raw.exit_exact(s,iN,r) = fE.exitflag;
            raw.exit_gaussian(s,iN,r) = fG.exitflag;

            % Store mobility estimates so tail behavior can be diagnosed exactly.
            raw.mu_hat_exact(s,iN,r,:) = reshape(fE.mu,1,1,1,K-1);
            raw.mu_hat_gaussian(s,iN,r,:) = reshape(fG.mu,1,1,1,K-1);
            raw.max_mu_exact(s,iN,r) = max(fE.mu);
            raw.max_mu_gaussian(s,iN,r) = max(fG.mu);

            % Keep a representative fitted transition matrix from the first
            % replication at each design point for qualitative diagnostics.
            if r==1
                raw.P_hat_exact{s,iN} = fE.P;
                raw.P_hat_gaussian{s,iN} = fG.P;
                raw.Y_example{s,iN} = Y;
            end

            completed(s,iN,r) = true;

            % Save after every replication so interruption loses at most one fit.
            save(checkpoint_file,'raw','completed','N_grid','T','R','delta_grid', ...
                'scenario_name','pi_true','mu_true','P_true','Delta','nstarts','h','-v7.3');

            fprintf(['scenario=%-15s N=%2d rep=%2d/%2d | cache=%6.3fs ' ...
                     'exact-opt=%6.3fs gauss-opt=%6.3fs | exit=(%d,%d)\n'], ...
                scenario_name{s},N,r,R,raw.cache_time(s,iN,r), ...
                raw.opt_time_exact(s,iN,r),raw.opt_time_gaussian(s,iN,r), ...
                fE.exitflag,fG.exitflag);
        end
    end
end

%% Summarize
template = struct( ...
    'scenario','', 'delta',NaN, 'N',NaN, ...
    'rmse_pi_exact',NaN, 'rmse_mu_exact',NaN, 'rmse_P_exact',NaN, ...
    'rmse_pi_gaussian',NaN, 'rmse_mu_gaussian',NaN, 'rmse_P_gaussian',NaN, ...
    'D_pi',NaN, 'D_mu',NaN, 'D_P',NaN, ...
    'cache_time_median',NaN, 'opt_time_exact_median',NaN, ...
    'end_time_exact_median',NaN, 'opt_time_gaussian_median',NaN, ...
    'end_over_gaussian',NaN, 'exact_success_rate',NaN, ...
    'gaussian_success_rate',NaN, ...
    'median_rel_mu_exact',NaN, 'q75_rel_mu_exact',NaN, 'q90_rel_mu_exact',NaN, ...
    'median_rel_mu_gaussian',NaN, 'q75_rel_mu_gaussian',NaN, 'q90_rel_mu_gaussian',NaN, ...
    'frac_maxmu_gt10_exact',NaN, 'frac_maxmu_gt10_gaussian',NaN, ...
    'P_hat_exact',[], 'P_hat_gaussian',[]);
results = repmat(template,nS,nN);
for s=1:nS
    for iN=1:nN
        z = template;
        z.scenario = scenario_name{s};
        z.delta = delta_grid(s);
        z.N = N_grid(iN);

        z.rmse_pi_exact = sqrt(mean(squeeze(raw.se_pi_exact(s,iN,:))));
        z.rmse_mu_exact = sqrt(mean(squeeze(raw.se_mu_exact(s,iN,:))));
        z.rmse_P_exact = sqrt(mean(squeeze(raw.se_P_exact(s,iN,:))));
        z.rmse_pi_gaussian = sqrt(mean(squeeze(raw.se_pi_gaussian(s,iN,:))));
        z.rmse_mu_gaussian = sqrt(mean(squeeze(raw.se_mu_gaussian(s,iN,:))));
        z.rmse_P_gaussian = sqrt(mean(squeeze(raw.se_P_gaussian(s,iN,:))));

        z.D_pi = sqrt(mean(squeeze(raw.d2_pi(s,iN,:))));
        z.D_mu = sqrt(mean(squeeze(raw.d2_mu(s,iN,:))));
        z.D_P = sqrt(mean(squeeze(raw.d2_P(s,iN,:))));

        z.cache_time_median = median(squeeze(raw.cache_time(s,iN,:)));
        z.opt_time_exact_median = median(squeeze(raw.opt_time_exact(s,iN,:)));
        z.end_time_exact_median = median(squeeze(raw.end_time_exact(s,iN,:)));
        z.opt_time_gaussian_median = median(squeeze(raw.opt_time_gaussian(s,iN,:)));
        z.end_over_gaussian = z.end_time_exact_median/z.opt_time_gaussian_median;

        z.exact_success_rate = mean(squeeze(raw.exit_exact(s,iN,:))>0);
        z.gaussian_success_rate = mean(squeeze(raw.exit_gaussian(s,iN,:))>0);

        % Robust mobility-tail diagnostics.
        muE = squeeze(raw.mu_hat_exact(s,iN,:,:));      % R-by-(K-1)
        muG = squeeze(raw.mu_hat_gaussian(s,iN,:,:));
        relE = sqrt(sum((muE-mu_true).^2,2))/norm(mu_true);
        relG = sqrt(sum((muG-mu_true).^2,2))/norm(mu_true);

        z.median_rel_mu_exact = median(relE,'omitnan');
        z.q75_rel_mu_exact = prctile(relE,75);
        z.q90_rel_mu_exact = prctile(relE,90);
        z.median_rel_mu_gaussian = median(relG,'omitnan');
        z.q75_rel_mu_gaussian = prctile(relG,75);
        z.q90_rel_mu_gaussian = prctile(relG,90);

        z.frac_maxmu_gt10_exact = mean(squeeze(raw.max_mu_exact(s,iN,:))>10,'omitnan');
        z.frac_maxmu_gt10_gaussian = mean(squeeze(raw.max_mu_gaussian(s,iN,:))>10,'omitnan');

        z.P_hat_exact = raw.P_hat_exact{s,iN};
        z.P_hat_gaussian = raw.P_hat_gaussian{s,iN};
        results(s,iN) = z;
    end
end

save(result_file,'results','raw','pi_true','mu_true','P_true','N_grid', ...
    'delta_grid','scenario_name','R','T','Delta','nstarts','h','-v7.3');

fprintf('\nExperiment 1 complete. Results written to %s\n',result_file);
% Plot after inspecting the pilot diagnostics.

%% Local helper
function raw = initialize_raw(nS,nN,R,K)
raw.se_pi_exact = nan(nS,nN,R);
raw.se_mu_exact = nan(nS,nN,R);
raw.se_P_exact = nan(nS,nN,R);
raw.se_pi_gaussian = nan(nS,nN,R);
raw.se_mu_gaussian = nan(nS,nN,R);
raw.se_P_gaussian = nan(nS,nN,R);
raw.d2_pi = nan(nS,nN,R);
raw.d2_mu = nan(nS,nN,R);
raw.d2_P = nan(nS,nN,R);
raw.cache_time = nan(nS,nN,R);
raw.opt_time_exact = nan(nS,nN,R);
raw.end_time_exact = nan(nS,nN,R);
raw.opt_time_gaussian = nan(nS,nN,R);
raw.exit_exact = nan(nS,nN,R);
raw.exit_gaussian = nan(nS,nN,R);
raw.mu_hat_exact = nan(nS,nN,R,K-1);
raw.mu_hat_gaussian = nan(nS,nN,R,K-1);
raw.max_mu_exact = nan(nS,nN,R);
raw.max_mu_gaussian = nan(nS,nN,R);
raw.P_hat_exact = cell(nS,nN);
raw.P_hat_gaussian = cell(nS,nN);
raw.Y_example = cell(nS,nN);
end

%% finalize_experiment1_R500_checkpoint.m
% Finalize Experiment 1 from an already completed checkpoint.
% This script DOES NOT rerun any Monte Carlo replications.
%
% Expected raw-array layout:
%   raw.field(scenario, N-index, replication)

clear; clc;

checkpoint_file = 'experiment1_checkpoint.mat';
result_file     = 'experiment1_results.mat';

if ~isfile(checkpoint_file)
    error('Cannot find %s in the current folder.', checkpoint_file);
end

S = load(checkpoint_file);
required = {'raw','completed','N_grid','T','R','delta_grid','scenario_name', ...
            'pi_true','mu_true','P_true','Delta','nstarts','h'};
for jj = 1:numel(required)
    if ~isfield(S,required{jj})
        error('Checkpoint is missing field %s.', required{jj});
    end
end

raw           = S.raw;
completed     = S.completed;
N_grid        = S.N_grid;
T             = S.T;
R             = S.R;
delta_grid    = S.delta_grid;
scenario_name = S.scenario_name;
pi_true       = S.pi_true;
mu_true       = S.mu_true;
P_true        = S.P_true;
Delta         = S.Delta;
nstarts       = S.nstarts;
h             = S.h;

nS = numel(delta_grid);
nN = numel(N_grid);

% Sanity checks: no simulation should be missing.
if ~isequal(size(completed), [nS nN R])
    error('Unexpected completed-array size. Expected [%d %d %d], got [%s].', ...
        nS,nN,R,num2str(size(completed)));
end
n_done = nnz(completed);
n_total = numel(completed);
fprintf('Checkpoint completion: %d/%d replications (%.2f%%).\n', ...
    n_done,n_total,100*n_done/n_total);
if n_done ~= n_total
    error('Checkpoint is not complete: %d replications are still missing.', n_total-n_done);
end

% Confirm the expected raw layout.
check_fields = {'se_pi_exact','se_mu_exact','se_P_exact', ...
                'se_pi_gaussian','se_mu_gaussian','se_P_gaussian', ...
                'd2_pi','d2_mu','d2_P','cache_time','opt_time_exact', ...
                'end_time_exact','opt_time_gaussian','exit_exact','exit_gaussian'};
for jj = 1:numel(check_fields)
    f = check_fields{jj};
    if ~isfield(raw,f)
        error('raw is missing field %s.',f);
    end
    sz = size(raw.(f));
    if ~isequal(sz,[nS nN R])
        error('raw.%s has size [%s], expected [%d %d %d].', ...
            f,num2str(sz),nS,nN,R);
    end
end

% Preallocate a homogeneous struct array explicitly. This avoids MATLAB's
% "Subscripted assignment between dissimilar structures" failure.
template = struct( ...
    'scenario','', ...
    'delta',NaN, ...
    'N',NaN, ...
    'rmse_pi_exact',NaN, ...
    'rmse_mu_exact',NaN, ...
    'rmse_P_exact',NaN, ...
    'rmse_pi_gaussian',NaN, ...
    'rmse_mu_gaussian',NaN, ...
    'rmse_P_gaussian',NaN, ...
    'D_pi',NaN, ...
    'D_mu',NaN, ...
    'D_P',NaN, ...
    'cache_time_median',NaN, ...
    'opt_time_exact_median',NaN, ...
    'end_time_exact_median',NaN, ...
    'opt_time_gaussian_median',NaN, ...
    'end_over_gaussian',NaN, ...
    'exact_success_rate',NaN, ...
    'gaussian_success_rate',NaN, ...
    'P_hat_exact',[], ...
    'P_hat_gaussian',[]);

results = repmat(template,nS,nN);

for s = 1:nS
    for iN = 1:nN
        z = template;
        z.scenario = scenario_name{s};
        z.delta = delta_grid(s);
        z.N = N_grid(iN);

        z.rmse_pi_exact = sqrt(mean(squeeze(raw.se_pi_exact(s,iN,:)),'omitnan'));
        z.rmse_mu_exact = sqrt(mean(squeeze(raw.se_mu_exact(s,iN,:)),'omitnan'));
        z.rmse_P_exact  = sqrt(mean(squeeze(raw.se_P_exact(s,iN,:)),'omitnan'));

        z.rmse_pi_gaussian = sqrt(mean(squeeze(raw.se_pi_gaussian(s,iN,:)),'omitnan'));
        z.rmse_mu_gaussian = sqrt(mean(squeeze(raw.se_mu_gaussian(s,iN,:)),'omitnan'));
        z.rmse_P_gaussian  = sqrt(mean(squeeze(raw.se_P_gaussian(s,iN,:)),'omitnan'));

        z.D_pi = sqrt(mean(squeeze(raw.d2_pi(s,iN,:)),'omitnan'));
        z.D_mu = sqrt(mean(squeeze(raw.d2_mu(s,iN,:)),'omitnan'));
        z.D_P  = sqrt(mean(squeeze(raw.d2_P(s,iN,:)),'omitnan'));

        z.cache_time_median       = median(squeeze(raw.cache_time(s,iN,:)),'omitnan');
        z.opt_time_exact_median   = median(squeeze(raw.opt_time_exact(s,iN,:)),'omitnan');
        z.end_time_exact_median   = median(squeeze(raw.end_time_exact(s,iN,:)),'omitnan');
        z.opt_time_gaussian_median= median(squeeze(raw.opt_time_gaussian(s,iN,:)),'omitnan');
        z.end_over_gaussian       = z.end_time_exact_median/z.opt_time_gaussian_median;

        z.exact_success_rate    = mean(squeeze(raw.exit_exact(s,iN,:))>0,'omitnan');
        z.gaussian_success_rate = mean(squeeze(raw.exit_gaussian(s,iN,:))>0,'omitnan');

        if isfield(raw,'P_hat_exact') && numel(raw.P_hat_exact) >= sub2ind([nS nN],s,iN)
            z.P_hat_exact = raw.P_hat_exact{s,iN};
        end
        if isfield(raw,'P_hat_gaussian') && numel(raw.P_hat_gaussian) >= sub2ind([nS nN],s,iN)
            z.P_hat_gaussian = raw.P_hat_gaussian{s,iN};
        end

        results(s,iN) = z;
    end
end

save(result_file,'results','raw','pi_true','mu_true','P_true','N_grid', ...
    'delta_grid','scenario_name','R','T','Delta','nstarts','h','-v7.3');

fprintf('\nExperiment 1 finalized successfully from checkpoint.\n');
fprintf('Saved %s with R=%d.\n',result_file,R);

% Compact summary to screen.
fprintf('\n%-17s %4s %11s %11s %11s %10s\n', ...
    'scenario','N','D_pi','D_mu','D_P','speedup');
for s=1:nS
    for iN=1:nN
        z=results(s,iN);
        fprintf('%-17s %4d %11.4g %11.4g %11.4g %10.2f\n', ...
            z.scenario,z.N,z.D_pi,z.D_mu,z.D_P,z.end_over_gaussian);
    end
end

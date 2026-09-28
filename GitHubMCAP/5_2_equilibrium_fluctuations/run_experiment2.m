% Public reproducibility script.
% Expected final result file: experiment2_results.mat
% After running, move/copy the result into the matching results/experiment*/ folder
% if the script saves it in the current MATLAB working directory.

%% experiment2_equilibrium_recovery.m
% Experiment 2: recovery from equilibrium fluctuations.
% Requires build_generator.m

clear; clc;
rng(20260907,'twister');

K = 4;
Delta = 1;
N = 1000;
pi_true = [0.12 0.28 0.38 0.22];
mu_base = [0.35 0.25 0.30];
scale_grid = [0.6 1.0 1.5];
M_grid = [50 100 250 500 1000 2000 5000];
R = 100;

nS = numel(scale_grid);
nM = numel(M_grid);
Mmax = max(M_grid);
y0 = round(N*pi_true);

if any(y0 ~= N*pi_true)
    error('N*pi_true must be integer-valued.');
end

Dpi = diag(pi_true);
Dpi_invsqrt = diag(1./sqrt(pi_true));

err_mu = nan(nS,nM,R);
err_P  = nan(nS,nM,R);
err_Q  = nan(nS,nM,R);
true_mu = nan(nS,K-1);
true_Q = nan(nS,K,K);
true_P = nan(nS,K,K);
eig_S2_true = nan(nS,K);
sensitivity = nan(nS,1);
theory_err_mu = nan(nS,1);
theory_err_P  = nan(nS,1);
theory_err_Q  = nan(nS,1);

for is = 1:nS
    scale = scale_grid(is);
    mu_true = scale*mu_base;
    Q = build_generator(pi_true,mu_true);
    P = expm(Delta*Q);

    true_mu(is,:) = mu_true;
    true_Q(is,:,:) = Q;
    true_P(is,:,:) = P;

    Sigma_true = N*(Dpi - P'*Dpi*P);
    C_true = (1/N)*Dpi_invsqrt*Sigma_true*Dpi_invsqrt;
    S2_true = eye(K) - C_true;
    S2_true = (S2_true + S2_true')/2;

    lam = sort(real(eig(S2_true)),'descend');
    eig_S2_true(is,:) = lam(:).';
    sensitivity(is) = 1/min(lam);

    [P_pop,Q_pop,mu_pop] = recover_from_covariance(Sigma_true,N,pi_true,Delta);
    theory_err_mu(is) = norm(mu_pop-mu_true);
    theory_err_P(is)  = norm(P_pop-P,'fro');
    theory_err_Q(is)  = norm(Q_pop-Q,'fro');

    fprintf('\nScale %.1f\n',scale);
    fprintf('  true mu = [%g %g %g]\n',mu_true);
    fprintf('  population mu error = %.3e\n',theory_err_mu(is));
    fprintf('  population P error  = %.3e\n',theory_err_P(is));
    fprintf('  population Q error  = %.3e\n',theory_err_Q(is));
    fprintf('  eig(S^2) ='); fprintf(' %.8g',lam); fprintf('\n');
    fprintf('  sensitivity = %.6g\n',sensitivity(is));

    row_cdf = cumsum(P,2);

    for r = 1:R
        Y = zeros(Mmax,K);

        for m = 1:Mmax
            y1 = zeros(1,K);
            for j = 1:K
                nj = y0(j);
                u = rand(nj,1);
                dest = 1 + sum(u > row_cdf(j,:),2);
                for k = 1:K
                    y1(k) = y1(k) + sum(dest==k);
                end
            end
            Y(m,:) = y1;
        end

        for iM = 1:nM
            M = M_grid(iM);
            Z = Y(1:M,:) - N*pi_true;
            Sigma_hat = (Z'*Z)/M;  % known equilibrium mean

            [P_hat,Q_hat,mu_hat] = recover_from_covariance(Sigma_hat,N,pi_true,Delta);

            err_mu(is,iM,r) = norm(mu_hat-mu_true);
            err_P(is,iM,r)  = norm(P_hat-P,'fro');
            err_Q(is,iM,r)  = norm(Q_hat-Q,'fro');
        end

        if mod(r,10)==0 || r==R
            fprintf('  completed %3d/%3d replications\n',r,R);
        end
    end
end

median_err_mu = median(err_mu,3,'omitnan');
median_err_P  = median(err_P,3,'omitnan');
median_err_Q  = median(err_Q,3,'omitnan');

q10_err_mu = quantile(err_mu,0.10,3); q90_err_mu = quantile(err_mu,0.90,3);
q10_err_P  = quantile(err_P,0.10,3);  q90_err_P  = quantile(err_P,0.90,3);
q10_err_Q  = quantile(err_Q,0.10,3);  q90_err_Q  = quantile(err_Q,0.90,3);

save('experiment2_results.mat', ...
    'K','Delta','N','pi_true','mu_base','scale_grid','M_grid','R','y0', ...
    'true_mu','true_Q','true_P','eig_S2_true','sensitivity', ...
    'theory_err_mu','theory_err_P','theory_err_Q', ...
    'err_mu','err_P','err_Q', ...
    'median_err_mu','median_err_P','median_err_Q', ...
    'q10_err_mu','q90_err_mu','q10_err_P','q90_err_P','q10_err_Q','q90_err_Q');

fprintf('\n============================================================\n');
fprintf('Experiment 2 complete. Saved experiment2_results.mat\n');
fprintf('============================================================\n');

for is = 1:nS
    fprintf('\nScale %.1f\n',scale_grid(is));
    fprintf('%7s | %12s %12s %12s\n','M','mu error','P error','Q error');
    for iM = 1:nM
        fprintf('%7d | %12.4g %12.4g %12.4g\n',M_grid(iM), ...
            median_err_mu(is,iM),median_err_P(is,iM),median_err_Q(is,iM));
    end
end

function [P_hat,Q_hat,mu_hat] = recover_from_covariance(Sigma_hat,N,pi_true,Delta)
    K = numel(pi_true);
    Dpi_sqrt = diag(sqrt(pi_true));
    Dpi_invsqrt = diag(1./sqrt(pi_true));

    C_hat = (1/N)*Dpi_invsqrt*Sigma_hat*Dpi_invsqrt;
    C_hat = (C_hat + C_hat')/2;
    S2_hat = eye(K)-C_hat;
    S2_hat = (S2_hat + S2_hat')/2;

    [U,L] = eig(S2_hat);
    lam = real(diag(L));
    lam = min(max(lam,1e-10),1);

    S_hat = U*diag(sqrt(lam))*U';
    H_hat = (1/(2*Delta))*U*diag(log(lam))*U';

    P_hat = Dpi_invsqrt*S_hat*Dpi_sqrt;
    Q_hat = Dpi_invsqrt*H_hat*Dpi_sqrt;

    mu_hat = zeros(1,K-1);
    for j = 1:K-1
        mu_hat(j) = H_hat(j,j+1);
    end
end

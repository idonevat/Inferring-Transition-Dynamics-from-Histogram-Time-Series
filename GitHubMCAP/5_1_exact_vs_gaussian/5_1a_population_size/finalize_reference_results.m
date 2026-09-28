%% finalize_experiment1_results.m
% Paper-ready summaries for Experiment 1.
% Requires experiment1_results.mat (or experiment1_results(1).mat).
%
% Produces:
%   experiment1_main_results.png
%   experiment1_agreement_runtime.png
%   experiment1_summary_table.csv
%   experiment1_runaway_table.csv
%
% Main figure:
%   rows = equilibrium / non-equilibrium
%   cols = RMSE(pi), median ||mu_hat-mu||_2, RMSE(P)
%
% Exact-vs-Gaussian agreement and runtime are summarized separately.

clear; clc;

if isfile('experiment1_results.mat')
    fname='experiment1_results.mat';
elseif isfile('experiment1_results(1).mat')
    fname='experiment1_results(1).mat';
else
    error('Cannot find experiment1_results.mat or experiment1_results(1).mat.');
end
S=load(fname);

raw=S.raw;
N_grid=S.N_grid(:).';
K=numel(S.pi_true);
nS=numel(S.delta_grid);
nN=numel(N_grid);
R=S.R;

% Error norms reconstructed from stored per-coordinate MSE quantities.
err_pi_E=sqrt(K*raw.se_pi_exact);
err_pi_G=sqrt(K*raw.se_pi_gaussian);
err_mu_E=sqrt((K-1)*raw.se_mu_exact);
err_mu_G=sqrt((K-1)*raw.se_mu_gaussian);
err_P_E=sqrt(K^2*raw.se_P_exact);
err_P_G=sqrt(K^2*raw.se_P_gaussian);

% Exact-vs-Gaussian disagreement norms.
d_pi=sqrt(raw.d2_pi);
d_mu=sqrt(raw.d2_mu);
d_P=sqrt(raw.d2_P);

% Robust runaway definition. There is a large empirical gap: ordinary
% errors are O(1), whereas runaway errors are >= O(1e5).
runaway_threshold=100;

% Verify threshold sensitivity.
thresholds=[10 100 1e3 1e5];
fprintf('\n=== Experiment 1 final summary ===\n');
fprintf('Runaway-count sensitivity:\n');
for q=1:numel(thresholds)
    ne=sum(err_mu_E(:)>thresholds(q));
    ng=sum(err_mu_G(:)>thresholds(q));
    fprintf(' threshold=%g: exact=%d, gaussian=%d\n',thresholds(q),ne,ng);
end

rows=[];
for s=1:nS
    for iN=1:nN
        z.scenario=s;
        z.delta=S.delta_grid(s);
        z.N=N_grid(iN);

        x=squeeze(err_pi_E(:,iN,s)); z.rmse_pi_exact=sqrt(mean(x.^2));
        x=squeeze(err_pi_G(:,iN,s)); z.rmse_pi_gaussian=sqrt(mean(x.^2));

        x=squeeze(err_mu_E(:,iN,s)); z.median_mu_exact=median(x);
        x=squeeze(err_mu_G(:,iN,s)); z.median_mu_gaussian=median(x);

        x=squeeze(err_P_E(:,iN,s)); z.rmse_P_exact=sqrt(mean(x.^2));
        x=squeeze(err_P_G(:,iN,s)); z.rmse_P_gaussian=sqrt(mean(x.^2));

        z.D_pi_median=median(squeeze(d_pi(:,iN,s)));
        z.D_mu_median=median(squeeze(d_mu(:,iN,s)));
        z.D_P_median=median(squeeze(d_P(:,iN,s)));

        z.runaway_exact=mean(squeeze(err_mu_E(:,iN,s))>runaway_threshold);
        z.runaway_gaussian=mean(squeeze(err_mu_G(:,iN,s))>runaway_threshold);

        z.cache_time_median=median(squeeze(raw.cache_time(:,iN,s)));
        z.exact_end_time_median=median(squeeze(raw.end_time_exact(:,iN,s)));
        z.gaussian_opt_time_median=median(squeeze(raw.opt_time_gaussian(:,iN,s)));

        rows=[rows; struct2table(z)]; %#ok<AGROW>
    end
end

writetable(rows,'experiment1_summary_table.csv');

runT=rows(:,{'scenario','delta','N','runaway_exact','runaway_gaussian'});
writetable(runT,'experiment1_runaway_table.csv');

disp(rows);

%% Main paper figure
figure('Color','w','Name','Experiment 1: exact versus Gaussian');
tl=tiledlayout(nS,3,'Padding','compact','TileSpacing','compact');

for s=1:nS
    q=rows(rows.scenario==s,:);

    nexttile;
    plot(q.N,q.rmse_pi_exact,'-o','LineWidth',1.4); hold on;
    plot(q.N,q.rmse_pi_gaussian,'--s','LineWidth',1.4);
    grid on; xlabel('N'); ylabel('RMSE(\pi)');
    if s==1, title('Composition'); end
    if s==1, legend('Exact','Gaussian','Location','best'); end

    nexttile;
    plot(q.N,q.median_mu_exact,'-o','LineWidth',1.4); hold on;
    plot(q.N,q.median_mu_gaussian,'--s','LineWidth',1.4);
    grid on; xlabel('N'); ylabel('Median ||\mu-hat-\mu||_2');
    if s==1, title('Mobility'); end

    nexttile;
    plot(q.N,q.rmse_P_exact,'-o','LineWidth',1.4); hold on;
    plot(q.N,q.rmse_P_gaussian,'--s','LineWidth',1.4);
    grid on; xlabel('N'); ylabel('RMSE(P)');
    if s==1, title('Transition matrix'); end
end
title(tl,'Experiment 1: exact likelihood versus Gaussian approximation');
exportgraphics(gcf,'experiment1_main_results.png','Resolution',250);

%% Agreement + runaway + runtime
figure('Color','w','Name','Experiment 1 diagnostics');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
for s=1:nS
    q=rows(rows.scenario==s,:);
    plot(q.N,q.D_pi_median,'-o','LineWidth',1.3); hold on;
end
grid on; xlabel('N'); ylabel('Median ||pi_E-pi_G||_2');
title('Exact-Gaussian composition agreement');
legend('Scenario 1','Scenario 2','Location','best');

nexttile;
for s=1:nS
    q=rows(rows.scenario==s,:);
    plot(q.N,q.D_mu_median,'-o','LineWidth',1.3); hold on;
end
grid on; xlabel('N'); ylabel('Median ||mu_E-mu_G||_2');
title('Exact-Gaussian mobility agreement');

nexttile;
for s=1:nS
    q=rows(rows.scenario==s,:);
    plot(q.N,100*q.runaway_exact,'-o','LineWidth',1.3); hold on;
    plot(q.N,100*q.runaway_gaussian,'--s','LineWidth',1.3);
end
grid on; xlabel('N'); ylabel('Runaway fits (%)');
title(sprintf('Practical fast-interface limit (threshold %g)',runaway_threshold));

nexttile;
q=rows(rows.scenario==1,:);
semilogy(q.N,q.exact_end_time_median,'-o','LineWidth',1.3); hold on;
semilogy(q.N,q.gaussian_opt_time_median,'--s','LineWidth',1.3);
grid on; xlabel('N'); ylabel('Median time (s)');
title('Computational cost');
legend('Exact end-to-end','Gaussian optimization','Location','best');

exportgraphics(gcf,'experiment1_agreement_runtime.png','Resolution',250);

fprintf('\nSaved paper-ready Experiment 1 summaries and figures.\n');
fprintf('Main mobility summary is median Euclidean error; RMSE(pi) and RMSE(P) remain conventional.\n');

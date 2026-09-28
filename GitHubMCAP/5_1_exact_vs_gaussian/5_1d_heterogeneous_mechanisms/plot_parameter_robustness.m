%% Experiment 1 robustness figure -- robust summary version
% Requires:
%   experiment1_parameter_robustness_results.mat
%
% Design:
%   Each panel shows medians and interquartile ranges (IQRs) across the
%   R=100 randomly generated mechanisms. Extreme flat-likelihood mobility
%   fits are deliberately NOT used to set the plot scale; their frequencies
%   should be reported separately in the manuscript text.
%
clear; close all; clc;

S = load('experiment1_parameter_robustness_results.mat');
N = S.N_grid(:)';
R = S.raw;
nN = numel(N);

% Convert all arrays to nN-by-Rpar (rows index population size).
relP_E = to_N_by_R(R.rel_P_exact,nN);
relP_G = to_N_by_R(R.rel_P_gaussian,nN);
dirP   = to_N_by_R(R.direct_P,nN);

relMu_E = to_N_by_R(R.rel_mu_exact,nN);
relMu_G = to_N_by_R(R.rel_mu_gaussian,nN);
dirMu   = to_N_by_R(R.direct_mu,nN);

% Robust summaries.
[p25_PE, med_PE, p75_PE] = robust_summary(relP_E);
[p25_PG, med_PG, p75_PG] = robust_summary(relP_G);
[p25_DP, med_DP, p75_DP] = robust_summary(dirP);

[p25_ME, med_ME, p75_ME] = robust_summary(relMu_E);
[p25_MG, med_MG, p75_MG] = robust_summary(relMu_G);
[p25_DM, med_DM, p75_DM] = robust_summary(dirMu);

% Slight horizontal offsets for paired exact / Gaussian summaries.
x = 1:nN;
dx = 0.075;

figure('Color','w','Position',[100 100 1120 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% (a) Error relative to truth for P
nexttile;
hold on;
plot_iqr(x-dx,med_PE,p25_PE,p75_PE,'o-');
plot_iqr(x+dx,med_PG,p25_PG,p75_PG,'s--');
hold off;
xticks(x); xticklabels(string(N));
xlabel('Population size N');
ylabel('Relative error in P');
title('(a) Estimation of the finite-time transition matrix');
legend({'Exact','Gaussian'},'Location','northeast');
grid on; box on;

%% (b) Direct Gaussian--exact discrepancy for P
nexttile;
plot_iqr(x,med_DP,p25_DP,p75_DP,'o-');
xticks(x); xticklabels(string(N));
xlabel('Population size N');
ylabel('Relative Gaussian--exact discrepancy in P');
title('(b) Approximation discrepancy for P');
grid on; box on;

%% (c) Error relative to truth for mobility
nexttile;
hold on;
plot_iqr(x-dx,med_ME,p25_ME,p75_ME,'o-');
plot_iqr(x+dx,med_MG,p25_MG,p75_MG,'s--');
hold off;
xticks(x); xticklabels(string(N));
xlabel('Population size N');
ylabel('Relative error in \mu');
title('(c) Estimation of transition intensities');
legend({'Exact','Gaussian'},'Location','northeast');
grid on; box on;

% Keep scale focused on the central 50% rather than pathological tails.
ylim([0, 1.25*max([p75_ME(:);p75_MG(:)])]);

%% (d) Direct Gaussian--exact discrepancy for mobility
nexttile;
plot_iqr(x,med_DM,p25_DM,p75_DM,'o-');
xticks(x); xticklabels(string(N));
xlabel('Population size N');
ylabel('Relative Gaussian--exact discrepancy in \mu');
title('(d) Approximation discrepancy for transition intensities');
grid on; box on;
ylim([0, 1.25*max(p75_DM)]);

% Global font settings.
set(findall(gcf,'-property','FontSize'),'FontSize',12);

exportgraphics(gcf,'experiment_1_robustness_summary.png','Resolution',300);

%% ---------- Local functions ----------
function Aout = to_N_by_R(A,nN)
    if size(A,1) == nN
        Aout = A;
    elseif size(A,2) == nN
        Aout = A.';
    else
        error('Unexpected array size %d-by-%d; neither dimension equals numel(N_grid)=%d.', ...
            size(A,1),size(A,2),nN);
    end
end

function [q25,med,q75] = robust_summary(A)
    q25 = prctile(A,25,2);
    med = median(A,2,'omitnan');
    q75 = prctile(A,75,2);
end

function plot_iqr(x,med,q25,q75,fmt)
    lo = med-q25;
    hi = q75-med;
    errorbar(x,med,lo,hi,fmt,'LineWidth',1.4,'MarkerSize',7, ...
        'CapSize',8);
end

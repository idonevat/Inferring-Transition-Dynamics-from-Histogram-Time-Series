% Public reproducibility script.
% Expected final result file: experiment3_results.mat
% After running, move/copy the result into the matching results/experiment*/ folder
% if the script saves it in the current MATLAB working directory.

%% experiment3_information_transition.m
% Experiment 3: Near-equilibrium information transition.
%
% Deterministic calculation of the Gaussian Fisher-information split
% for a one-dimensional mobility-scale direction
%
%     mu(beta) = exp(beta) * mu_true,
%
% evaluated at beta=0.  For
%
%     rho_N = pi + c N^{-alpha} h,
%
% the theorem predicts
%
%     I_mean  ~ N^(1-2 alpha),
%     I_fluct ~ O(1),
%
% while at exact equilibrium I_mean = 0.
%
% Outputs:
%   experiment3_results.mat
%   experiment3_information_transition.png
%   experiment3_slopes.csv

clear; clc;

%% Model settings
pi_true = [0.12 0.28 0.38 0.22];
mu_true = [0.35 0.25 0.30];
Delta = 1;

% Same perturbation direction used in Experiment 1.
h = [0.12 -0.06 -0.04 -0.02];
h = h - mean(h); % numerical safeguard: sum(h)=0
c = 1;

alpha_grid = 0.20:0.10:0.80;
N_grid = [50 100 200 500 1000 2000 5000 10000 20000];

K = numel(pi_true);
d = K-1;

% Central finite difference in log mobility scale beta.
eps_beta = 1e-5;

%% Allocate
nA = numel(alpha_grid);
nN = numel(N_grid);

Imean = nan(nA,nN);
Ifluct = nan(nA,nN);
Itotal = nan(nA,nN);
rho_store = nan(nA,nN,K);

% Exact equilibrium benchmark.
Imean_eq = zeros(1,nN);
Ifluct_eq = nan(1,nN);
Itotal_eq = nan(1,nN);

%% Helper: reduced Gaussian moments for given rho and beta
moment_fun = @(rho,beta) reduced_moments(rho, exp(beta)*mu_true, pi_true, Delta);

%% Local sequences
for ia = 1:nA
    alpha = alpha_grid(ia);

    for iN = 1:nN
        N = N_grid(iN);

        rho = pi_true + c*N^(-alpha)*h;
        if any(rho <= 0)
            error('rho_N left the simplex at alpha=%g, N=%g.',alpha,N);
        end
        rho = rho/sum(rho);
        rho_store(ia,iN,:) = rho;

        [m0,V0] = moment_fun(rho,0);
        [mp,Vp] = moment_fun(rho,+eps_beta);
        [mm,Vm] = moment_fun(rho,-eps_beta);

        dm = (mp-mm)/(2*eps_beta);
        dV = (Vp-Vm)/(2*eps_beta);

        V0 = (V0+V0')/2;
        dV = (dV+dV')/2;

        % Stable solves rather than explicit inverse.
        z = V0\dm;
        A = V0\dV;

        Imean(ia,iN) = N*(dm'*z);
        Ifluct(ia,iN) = 0.5*trace(A*A);
        Itotal(ia,iN) = Imean(ia,iN)+Ifluct(ia,iN);
    end
end

%% Exact equilibrium
rho = pi_true;
for iN = 1:nN
    N = N_grid(iN);

    [~,V0] = moment_fun(rho,0);
    [~,Vp] = moment_fun(rho,+eps_beta);
    [~,Vm] = moment_fun(rho,-eps_beta);

    dV = (Vp-Vm)/(2*eps_beta);
    V0 = (V0+V0')/2;
    dV = (dV+dV')/2;

    A = V0\dV;

    % The mean derivative is exactly zero at equilibrium for mobility.
    Imean_eq(iN) = 0;
    Ifluct_eq(iN) = 0.5*trace(A*A);
    Itotal_eq(iN) = Ifluct_eq(iN);
end

%% Estimate asymptotic slopes from the largest five N values
tail = max(1,nN-4):nN;
slope_mean = nan(nA,1);
slope_fluct = nan(nA,1);

for ia=1:nA
    p = polyfit(log(N_grid(tail)),log(Imean(ia,tail)),1);
    slope_mean(ia)=p(1);

    p = polyfit(log(N_grid(tail)),log(Ifluct(ia,tail)),1);
    slope_fluct(ia)=p(1);
end

expected_slope = 1-2*alpha_grid(:);

fprintf('\n=== Experiment 3: near-equilibrium information transition ===\n');
for ia=1:nA
    fprintf('alpha=%.2f: expected mean slope=%+.3f, estimated=%+.3f; fluctuation slope=%+.3f\n', ...
        alpha_grid(ia),expected_slope(ia),slope_mean(ia),slope_fluct(ia));
end

fprintf('\nExact equilibrium fluctuation information:\n');
fprintf('  min = %.8g, max = %.8g, relative range = %.3e\n', ...
    min(Ifluct_eq),max(Ifluct_eq), ...
    (max(Ifluct_eq)-min(Ifluct_eq))/mean(Ifluct_eq));

%% Save numerical results
slopes = table(alpha_grid(:),expected_slope,slope_mean,slope_fluct, ...
    'VariableNames',{'alpha','expected_mean_slope','estimated_mean_slope','estimated_fluctuation_slope'});
writetable(slopes,'experiment3_slopes.csv');

save('experiment3_results.mat', ...
    'pi_true','mu_true','Delta','h','c','alpha_grid','N_grid', ...
    'Imean','Ifluct','Itotal','Imean_eq','Ifluct_eq','Itotal_eq', ...
    'rho_store','slope_mean','slope_fluct','expected_slope','eps_beta');

%% Main-paper figure: information-channel transition versus N
% For each alpha, show the fraction of total Fisher information carried by
% the mean:
%
%     share_mean = I_mean / (I_mean + I_fluct).
%
% The seven alpha values are the same values used throughout the experiment.
% The critical case alpha = 1/2 is emphasized in the figure.
%
%   alpha < 1/2 : mean-information share increases with N
%   alpha = 1/2 : mean-information share remains O(1)
%   alpha > 1/2 : mean-information share decreases toward zero
%
% At exact equilibrium, mean information is identically zero.

mean_share = Imean ./ Itotal;

fig = figure( ...
    'Color','w', ...
    'Name','Experiment 3: Information-channel transition', ...
    'Units','inches', ...
    'Position',[1 1 5.1 3.55]);

ax = axes(fig);
hold(ax,'on');

% Plot all alpha values.  Use restrained line styling; alpha=0.5 is
% emphasized because it is the theoretical transition point.
for ia = 1:nA
    if abs(alpha_grid(ia)-0.5) < 1e-12
        semilogx(ax,N_grid,mean_share(ia,:),'-s', ...
            'LineWidth',2.2,'MarkerSize',5.2, ...
            'MarkerIndices',1:2:numel(N_grid), ...
            'DisplayName',sprintf('\\alpha = %.1f  (critical)',alpha_grid(ia)));
    else
        semilogx(ax,N_grid,mean_share(ia,:),'-o', ...
            'LineWidth',1.25,'MarkerSize',3.8, ...
            'MarkerIndices',1:2:numel(N_grid), ...
            'DisplayName',sprintf('\\alpha = %.1f',alpha_grid(ia)));
    end
end

% Exact-equilibrium reference: deliberately unobtrusive.
semilogx(ax,N_grid,zeros(size(N_grid)),'--', ...
    'LineWidth',1.0, ...
    'DisplayName','exact equilibrium');

hold(ax,'off');

xlabel(ax,'Population size, N');
ylabel(ax,'Fraction of information carried by the mean');
xlim(ax,[min(N_grid) max(N_grid)]);
ylim(ax,[-0.025 1]);
yticks(ax,0:0.2:1);

grid(ax,'on');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
box(ax,'on');

legend(ax,'Location','eastoutside','Box','off');

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',9.5, ...
    'LineWidth',0.8, ...
    'TickDir','out');

% Helpful annotations identifying the two asymptotic directions.
text(ax,0.985,0.94,'mean dominated', ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontAngle','italic', ...
    'FontSize',8.5);

text(ax,0.985,0.055,'fluctuation dominated', ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', ...
    'FontAngle','italic', ...
    'FontSize',8.5);

exportgraphics(fig,'experiment3_information_channels_v3.pdf','ContentType','vector');
exportgraphics(fig,'experiment3_information_channels_v3.png','Resolution',300);
savefig(fig,'experiment3_information_channels_v3.fig');

fprintf('\nSaved main-paper figure:\n');
fprintf('  experiment3_information_channels_v3.pdf\n');
fprintf('  experiment3_information_channels_v3.png\n');
fprintf('  experiment3_information_channels_v3.fig\n');

%% Numerical summary for all plotted alpha values
fprintf('\nSlope verification over the full alpha grid:\n');
fprintf(' alpha | theoretical | estimated mean | estimated fluctuation\n');
fprintf('-------------------------------------------------------------\n');
for ia = 1:nA
    expected = 1 - 2*alpha_grid(ia);
    fprintf(' %.2f  |   %+7.3f   |    %+7.3f    |       %+7.3f\n', ...
        alpha_grid(ia),expected,slope_mean(ia),slope_fluct(ia));
end

fprintf('\nExact equilibrium fluctuation information:\n');
fprintf('  min = %.8f, max = %.8f, relative range = %.3e\n', ...
    min(Ifluct_eq),max(Ifluct_eq), ...
    (max(Ifluct_eq)-min(Ifluct_eq))/mean(Ifluct_eq));

fprintf('\nMean-information share I_mean / I_total:\n');
fprintf('%9s','N');
for ia = 1:nA
    fprintf(' | a=%.2f',alpha_grid(ia));
end
fprintf(' | equilibrium\n');

for iN = 1:nN
    fprintf('%9d',N_grid(iN));
    for ia = 1:nA
        fprintf(' | %6.4f',mean_share(ia,iN));
    end
    fprintf(' | %11.6f\n',0);
end

%% ------------------------------------------------------------------------
function [m,V] = reduced_moments(rho,mu,pi_vec,Delta)
% Reduced conditional mean proportion and covariance-per-individual.
%
% For one aggregate transition with source composition rho:
%   m = (rho P)_-
%   V = sum_j rho_j [diag(q_j)-q_j q_j']_-
% where q_j is row j of P.

K = numel(pi_vec);
Q = build_generator(pi_vec,mu);
P = expm(Delta*Q);

pnext = rho*P;
m = pnext(1:K-1).';

Vfull = zeros(K,K);
for j=1:K
    q = P(j,:).';
    Vfull = Vfull + rho(j)*(diag(q)-q*q.');
end
V = Vfull(1:K-1,1:K-1);
V = (V+V')/2;
end

%% Experiment 1: exact likelihood vs Gaussian approximation
% 4-panel publication figure using the saved R=200 results.
%
% Saved array layout in MATLAB:
%   raw.mu_hat_exact      : scenario x N x replication x mobility
%   raw.mu_hat_gaussian   : scenario x N x replication x mobility
%   raw.se_pi_*           : scenario x N x replication
%   raw.se_P_*            : scenario x N x replication
%
% The figure pools the two source-composition regimes at each N.
%
% Panels:
%   (a) median Euclidean error in pi
%   (b) median Frobenius error in P
%   (c) median relative Euclidean error in mu
%   (d) percentage of fits with max_j mu_hat_j > 10

clear; clc; close all;

%% Load results
dataFile = 'experiment1_reference_m93_R200_results.mat';
S = load(dataFile);

Ngrid = S.N_grid(:).';
mu0   = S.mu_true(:).';        % 1 x (K-1)
nmu   = norm(mu0);

nScenario = size(S.raw.mu_hat_exact,1);
nN        = size(S.raw.mu_hat_exact,2);
R         = size(S.raw.mu_hat_exact,3);
nMu       = size(S.raw.mu_hat_exact,4);

assert(nN == numel(Ngrid), ...
    'N_grid does not match the second dimension of mu_hat_exact.');
assert(nMu == numel(mu0), ...
    'mu_true does not match the fourth dimension of mu_hat_exact.');

fprintf('Detected layout: %d scenarios x %d N-values x %d replications x %d mobility parameters.\n', ...
    nScenario,nN,R,nMu);

%% Storage
medPiExact  = nan(1,nN);
medPiGauss  = nan(1,nN);

medPExact   = nan(1,nN);
medPGauss   = nan(1,nN);

medMuExact  = nan(1,nN);
medMuGauss  = nan(1,nN);

pathExact   = nan(1,nN);
pathGauss   = nan(1,nN);

%% Calculate pooled summaries
for i = 1:nN

    % -------------------------------------------------------------
    % pi error
    % raw.se_pi_* stores squared Euclidean estimation error.
    % Pool scenarios and replications before taking the median.
    % -------------------------------------------------------------
    ePiE = sqrt(reshape(S.raw.se_pi_exact(:,i,:),[],1));
    ePiG = sqrt(reshape(S.raw.se_pi_gaussian(:,i,:),[],1));

    medPiExact(i) = median(ePiE,'omitnan');
    medPiGauss(i) = median(ePiG,'omitnan');

    % -------------------------------------------------------------
    % P error
    % raw.se_P_* stores squared Frobenius estimation error.
    % -------------------------------------------------------------
    ePE = sqrt(reshape(S.raw.se_P_exact(:,i,:),[],1));
    ePG = sqrt(reshape(S.raw.se_P_gaussian(:,i,:),[],1));

    medPExact(i) = median(ePE,'omitnan');
    medPGauss(i) = median(ePG,'omitnan');

    % -------------------------------------------------------------
    % mu relative error and pathological-fit indicator
    %
    % For fixed scenario s and N index i:
    % squeeze(raw.mu_hat_*(s,i,:,:)) is R x (K-1).
    % -------------------------------------------------------------
    eMuE = [];
    eMuG = [];
    badE = [];
    badG = [];

    for s = 1:nScenario

        AE = squeeze(S.raw.mu_hat_exact(s,i,:,:));       % R x (K-1)
        AG = squeeze(S.raw.mu_hat_gaussian(s,i,:,:));    % R x (K-1)

        % Defensive orientation check for the R=1 edge case.
        if size(AE,2) ~= numel(mu0)
            AE = AE.';
        end
        if size(AG,2) ~= numel(mu0)
            AG = AG.';
        end

        % Relative Euclidean error in mobility vector
        errE = vecnorm(AE - mu0,2,2) / nmu;
        errG = vecnorm(AG - mu0,2,2) / nmu;

        eMuE = [eMuE; errE]; %#ok<AGROW>
        eMuG = [eMuG; errG]; %#ok<AGROW>

        % Runaway / fast-tail fit
        badE = [badE; max(AE,[],2) > 10]; %#ok<AGROW>
        badG = [badG; max(AG,[],2) > 10]; %#ok<AGROW>
    end

    medMuExact(i) = median(eMuE,'omitnan');
    medMuGauss(i) = median(eMuG,'omitnan');

    pathExact(i) = 100*mean(badE,'omitnan');
    pathGauss(i) = 100*mean(badG,'omitnan');
end

%% Print numerical values used in the figure
fprintf('\nPooled summaries over both regimes:\n');
fprintf('N        :'); fprintf(' %8g',Ngrid); fprintf('\n');

fprintf('pi exact :'); fprintf(' %8.4f',medPiExact); fprintf('\n');
fprintf('pi gauss :'); fprintf(' %8.4f',medPiGauss); fprintf('\n');

fprintf('P exact  :'); fprintf(' %8.4f',medPExact); fprintf('\n');
fprintf('P gauss  :'); fprintf(' %8.4f',medPGauss); fprintf('\n');

fprintf('mu exact :'); fprintf(' %8.4f',medMuExact); fprintf('\n');
fprintf('mu gauss :'); fprintf(' %8.4f',medMuGauss); fprintf('\n');

fprintf('bad exact:'); fprintf(' %8.2f',pathExact); fprintf('\n');
fprintf('bad gauss:'); fprintf(' %8.2f',pathGauss); fprintf('\n');

%% Plot
fig = figure('Color','w','Units','centimeters','Position',[2 2 23 19]);

tl = tiledlayout(fig,2,2,...
    'TileSpacing','compact',...
    'Padding','compact');

lw = 1.5;
ms = 5.5;

% ----- (a) equilibrium composition -------------------------------
ax1 = nexttile(tl,1);
hold(ax1,'on');
plot(ax1,Ngrid,medPiExact,'-o','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Exact likelihood');
plot(ax1,Ngrid,medPiGauss,'-s','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Gaussian approximation');
grid(ax1,'on');
box(ax1,'on');
xlabel(ax1,'Population size, $N$','Interpreter','latex');
ylabel(ax1,'Median error in $\pi$','Interpreter','latex');
title(ax1,'(a) Equilibrium composition','FontWeight','normal');
xlim(ax1,[min(Ngrid) max(Ngrid)]);

% ----- (b) transition matrix -------------------------------------
ax2 = nexttile(tl,2);
hold(ax2,'on');
plot(ax2,Ngrid,medPExact,'-o','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Exact likelihood');
plot(ax2,Ngrid,medPGauss,'-s','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Gaussian approximation');
grid(ax2,'on');
box(ax2,'on');
xlabel(ax2,'Population size, $N$','Interpreter','latex');
ylabel(ax2,'Median error in $P$','Interpreter','latex');
title(ax2,'(b) Transition matrix','FontWeight','normal');
xlim(ax2,[min(Ngrid) max(Ngrid)]);

% ----- (c) mobility ----------------------------------------------
ax3 = nexttile(tl,3);
hold(ax3,'on');
plot(ax3,Ngrid,medMuExact,'-o','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Exact likelihood');
plot(ax3,Ngrid,medMuGauss,'-s','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Gaussian approximation');
grid(ax3,'on');
box(ax3,'on');
xlabel(ax3,'Population size, $N$','Interpreter','latex');
ylabel(ax3,'Median relative error in $\mu$','Interpreter','latex');
title(ax3,'(c) Interface mobility','FontWeight','normal');
xlim(ax3,[min(Ngrid) max(Ngrid)]);

% ----- (d) pathological mobility estimates -----------------------

ax4 = nexttile(tl,4);
hold(ax4,'on');
plot(ax4,Ngrid,pathExact,'-o','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Exact likelihood');
plot(ax4,Ngrid,pathGauss,'-s','LineWidth',lw,'MarkerSize',ms,...
    'DisplayName','Gaussian approximation');
grid(ax4,'on');
box(ax4,'on');
xlabel(ax4,'Population size, $N$','Interpreter','latex');
ylabel(ax4,'Fits with $\max_j\hat{\mu}_j>10$ (\%)','Interpreter','latex');
title(ax4,'(d) Fast-tail estimates','FontWeight','normal');
xlim(ax4,[min(Ngrid) max(Ngrid)]);
ylim(ax4,[0 max([pathExact pathGauss])*1.12]);

% Common legend
lgd = legend(ax1,'Location','northoutside','Orientation','horizontal');
lgd.Layout.Tile = 'north';

% Consistent typography
allAxes = [ax1 ax2 ax3 ax4];
set(allAxes,...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'LineWidth',1,...
    'TickDir','out');

%% Export
exportgraphics(fig,'experiment1_reference_4panel.pdf',...
    'ContentType','vector');

exportgraphics(fig,'experiment1_reference_4panel.png',...
    'Resolution',400);

fprintf('\nSaved:\n');
fprintf('  experiment1_reference_4panel.pdf\n');
fprintf('  experiment1_reference_4panel.png\n');

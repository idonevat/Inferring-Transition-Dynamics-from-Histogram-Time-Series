%% plot_nonlocal_tree_robustness.m
% Plot the full nonlocal-tree robustness experiment.
%
% Expected input:
%   nonlocal_tree_robustness_full.mat
%
% Figure:
%   (a) equilibrium composition error E_pi
%   (b) mobility relative error E_mu
%   (c) transition-matrix error E_P
%
% Each panel shows the median across Monte Carlo replications, with
% 25th--75th percentile (IQR) error bars, for the nearest-neighbour and nonlocal models.
%
% The plotting order matches the numerical figures in the paper:
% E_pi -> E_mu -> E_P.

clear; clc; close all;

%% ---------------- SETTINGS ----------------------------------------------
results_file = 'nonlocal_tree_robustness_full.mat';

% Output names
fig_file = 'nonlocal_tree_robustness_figure.fig';
pdf_file = 'nonlocal_tree_robustness_figure.pdf';
png_file = 'nonlocal_tree_robustness_figure.png';

%% ---------------- LOAD RESULTS ------------------------------------------
S = load(results_file);

if isfield(S,'results')
    R = S.results;
else
    R = S;
end

T_grid = R.T_grid(:);

% Dimensions in saved file:
%   err_*(replication, T index, graph)
% graph 1 = nearest-neighbour
% graph 2 = nonlocal
err_pi = R.err_pi;
err_mu = R.err_mu;
err_P  = R.err_P;

graph_labels = {'Nearest-neighbour','Nonlocal'};

%% ---------------- SUMMARIES ---------------------------------------------
% median and IQR across replications (dimension 1)
med_pi = squeeze(median(err_pi,1,'omitnan'));
med_mu = squeeze(median(err_mu,1,'omitnan'));
med_P  = squeeze(median(err_P ,1,'omitnan'));

q25_pi = squeeze(prctile(err_pi,25,1));
q75_pi = squeeze(prctile(err_pi,75,1));

q25_mu = squeeze(prctile(err_mu,25,1));
q75_mu = squeeze(prctile(err_mu,75,1));

q25_P = squeeze(prctile(err_P,25,1));
q75_P = squeeze(prctile(err_P,75,1));

%% ---------------- FIGURE ------------------------------------------------
fig = figure('Color','w','Position',[100 100 1180 350]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

lw = 1.6;
ms = 6;

% MATLAB default axes color order; use the first two entries consistently.
co = get(groot,'defaultAxesColorOrder');
c1 = co(1,:);
c2 = co(2,:);

% Helper for an IQR band. Handle is hidden so the legend contains curves only.
band = @(x,lo,hi,c) fill([x;flipud(x)], [lo;flipud(hi)], c, ...
    'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');

% -------------------------------------------------------------------------
% (a) Equilibrium composition error
% -------------------------------------------------------------------------
ax1 = nexttile; hold(ax1,'on');
band(T_grid,q25_pi(:,1),q75_pi(:,1),c1);
band(T_grid,q25_pi(:,2),q75_pi(:,2),c2);
h1 = plot(T_grid,med_pi(:,1),'-o','LineWidth',lw,'MarkerSize',ms,'Color',c1);
h2 = plot(T_grid,med_pi(:,2),'--s','LineWidth',lw,'MarkerSize',ms,'Color',c2);
xlabel('Observation length, $T$','Interpreter','latex');
ylabel('Median error in $\pi$','Interpreter','latex');
grid on; box on;
set(gca,'XScale','log','XTick',T_grid,'XTickLabel',string(T_grid));

% -------------------------------------------------------------------------
% (b) Mobility relative error
% -------------------------------------------------------------------------
ax2 = nexttile; hold(ax2,'on');
band(T_grid,q25_mu(:,1),q75_mu(:,1),c1);
band(T_grid,q25_mu(:,2),q75_mu(:,2),c2);
plot(T_grid,med_mu(:,1),'-o','LineWidth',lw,'MarkerSize',ms,'Color',c1);
plot(T_grid,med_mu(:,2),'--s','LineWidth',lw,'MarkerSize',ms,'Color',c2);
xlabel('Observation length, $T$','Interpreter','latex');
ylabel('Median relative error in $\mu$','Interpreter','latex');
grid on; box on;
set(gca,'XScale','log','XTick',T_grid,'XTickLabel',string(T_grid));

% -------------------------------------------------------------------------
% (c) Transition-matrix error
% -------------------------------------------------------------------------
ax3 = nexttile; hold(ax3,'on');
band(T_grid,q25_P(:,1),q75_P(:,1),c1);
band(T_grid,q25_P(:,2),q75_P(:,2),c2);
plot(T_grid,med_P(:,1),'-o','LineWidth',lw,'MarkerSize',ms,'Color',c1);
plot(T_grid,med_P(:,2),'--s','LineWidth',lw,'MarkerSize',ms,'Color',c2);
xlabel('Observation length, $T$','Interpreter','latex');
ylabel('Median error in $P$','Interpreter','latex');
grid on; box on;
set(gca,'XScale','log','XTick',T_grid,'XTickLabel',string(T_grid));

%% ---------------- SHARED LEGEND -----------------------------------------
lgd = legend(ax1,[h1 h2],graph_labels, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
try
    lgd.Layout.Tile = 'north';
catch
end

%% ---------------- TYPOGRAPHY --------------------------------------------
set([ax1 ax2 ax3], ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'LineWidth',1, ...
    'TickDir','out', ...
    'TickLabelInterpreter','latex');

% Match the manuscript figure: major grid only.
for ax = [ax1 ax2 ax3]
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    xlim(ax,[min(T_grid) max(T_grid)]);
end

% Panel labels, positioned as in the T-sensitivity manuscript figure.
text(ax1,0.90,0.96,'(a)','Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',14);
text(ax2,0.90,0.96,'(b)','Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',14);
text(ax3,0.90,0.96,'(c)','Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontName','Times New Roman','FontSize',14);

% The shaded regions are the empirical 25th--75th percentile ranges across
% the R Monte Carlo replications. They represent dataset-to-dataset
% variability, not uncertainty in the median.

%% ---------------- PRINT NUMERICAL SUMMARY -------------------------------
fprintf('\nNonlocal-tree robustness: medians [IQR]\n');
fprintf('---------------------------------------------------------------\n');

for it = 1:numel(T_grid)
    fprintf('T = %d\n',T_grid(it));

    fprintf('  E_pi  Nearest-neighbour %.4f [%.4f, %.4f]   Nonlocal %.4f [%.4f, %.4f]\n', ...
        med_pi(it,1),q25_pi(it,1),q75_pi(it,1), ...
        med_pi(it,2),q25_pi(it,2),q75_pi(it,2));

    fprintf('  E_mu  Nearest-neighbour %.4f [%.4f, %.4f]   Nonlocal %.4f [%.4f, %.4f]\n', ...
        med_mu(it,1),q25_mu(it,1),q75_mu(it,1), ...
        med_mu(it,2),q25_mu(it,2),q75_mu(it,2));

    fprintf('  E_P   Nearest-neighbour %.4f [%.4f, %.4f]   Nonlocal %.4f [%.4f, %.4f]\n', ...
        med_P(it,1),q25_P(it,1),q75_P(it,1), ...
        med_P(it,2),q25_P(it,2),q75_P(it,2));
end

%% ---------------- SAVE --------------------------------------------------
savefig(fig,fig_file);

exportgraphics(fig,pdf_file,'ContentType','vector');
exportgraphics(fig,png_file,'Resolution',300);

fprintf('\nSaved:\n');
fprintf('  %s\n',fig_file);
fprintf('  %s\n',pdf_file);
fprintf('  %s\n',png_file);

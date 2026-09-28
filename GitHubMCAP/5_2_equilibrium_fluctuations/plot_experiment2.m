%% plot_experiment2_paper_figure_v2.m
% Simplified publication figure for Experiment 2
%
% Reads:
%   experiment2_results.mat
%
% Requires:
%   build_generator.m
%
% Panels:
%   (a) median relative recovery error in P versus M
%   (b) median relative recovery error in mu versus M
%   (c) conditioning diagnostic lambda_min(S^2) over a dense mobility grid
%
% No additional Monte Carlo simulation is required.

clear; clc;

%% Load completed Experiment 2
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
fname = fullfile(repoRoot,'results','experiment2','experiment2_results.mat');
if ~isfile(fname)
    error('Cannot find %s in the current folder.',fname);
end
S = load(fname);

M_grid     = S.M_grid(:).';
scale_grid = S.scale_grid(:).';
pi_true    = S.pi_true(:).';
mu_base    = S.mu_base(:).';
Delta      = S.Delta;

nS = numel(scale_grid);

%% Relative errors from the stored Monte Carlo results
% err_P and err_mu are absolute norm errors.  Normalize replication by
% replication using the corresponding fixed true norm for each scale.

rel_P  = nan(size(S.err_P));
rel_mu = nan(size(S.err_mu));

for is = 1:nS
    Ptrue  = squeeze(S.true_P(is,:,:));
    mutrue = S.true_mu(is,:);

    rel_P(is,:,:)  = S.err_P(is,:,:)  ./ norm(Ptrue,'fro');
    rel_mu(is,:,:) = S.err_mu(is,:,:) ./ norm(mutrue);
end

med_rel_P  = median(rel_P,3,'omitnan');
med_rel_mu = median(rel_mu,3,'omitnan');

%% Dense deterministic conditioning curve
% Under reversibility:
%   S = D_pi^(1/2) P D_pi^(-1/2) = exp(Delta H)
% and S^2 is positive definite.  Recovery through sqrt(S^2) and log(S^2)
% becomes more sensitive as lambda_min(S^2) approaches zero.

c_grid = linspace(0.2,2.0,181);
lambda_min = nan(size(c_grid));

Dpi_sqrt    = diag(sqrt(pi_true));
Dpi_invsqrt = diag(1./sqrt(pi_true));

for ic = 1:numel(c_grid)
    mu = c_grid(ic)*mu_base;
    Q  = build_generator(pi_true,mu);
    P  = expm(Delta*Q);

    Smat = Dpi_sqrt * P * Dpi_invsqrt;
    S2   = Smat*Smat;
    S2   = (S2+S2')/2;

    lambda_min(ic) = min(eig(S2));
end

% Exact locations of the three simulated mobility scales.
lambda_sim = nan(size(scale_grid));
for is = 1:nS
    mu = scale_grid(is)*mu_base;
    Q  = build_generator(pi_true,mu);
    P  = expm(Delta*Q);
    Smat = Dpi_sqrt * P * Dpi_invsqrt;
    S2 = (Smat*Smat + (Smat*Smat)')/2;
    lambda_sim(is) = min(eig(S2));
end

%% Print values used in the figure
fprintf('\nMedian relative recovery error in P\n');
fprintf('%7s', 'M');
for is=1:nS, fprintf(' | c=%.1f',scale_grid(is)); end
fprintf('\n');
for iM=1:numel(M_grid)
    fprintf('%7d',M_grid(iM));
    for is=1:nS
        fprintf(' | %6.4f',med_rel_P(is,iM));
    end
    fprintf('\n');
end

fprintf('\nMedian relative recovery error in mu\n');
fprintf('%7s', 'M');
for is=1:nS, fprintf(' | c=%.1f',scale_grid(is)); end
fprintf('\n');
for iM=1:numel(M_grid)
    fprintf('%7d',M_grid(iM));
    for is=1:nS
        fprintf(' | %6.4f',med_rel_mu(is,iM));
    end
    fprintf('\n');
end

fprintf('\nConditioning at simulated scales\n');
for is=1:nS
    fprintf('c = %.1f: lambda_min(S^2) = %.6f\n', ...
        scale_grid(is),lambda_sim(is));
end

%% Figure
fig = figure( ...
    'Units','inches', ...
    'Position',[1 1 7.2 3.0], ...
    'Color','w');

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

labels = arrayfun(@(c) sprintf('c = %.1f',c), ...
    scale_grid,'UniformOutput',false);

%% (a) P
ax1 = nexttile(tl,1);
hold(ax1,'on');
for is = 1:nS
    plot(ax1,M_grid,med_rel_P(is,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',4.5, ...
        'DisplayName',labels{is});
end
hold(ax1,'off');

set(ax1,'XScale','log','YScale','log');
xlabel(ax1,'Number of transitions, M');
ylabel(ax1,'Median relative error');
title(ax1,'(a) Recovery of P');
grid(ax1,'on');
box(ax1,'on');
legend(ax1,'Location','southwest','Box','off');

%% (b) mu
ax2 = nexttile(tl,2);
hold(ax2,'on');
for is = 1:nS
    plot(ax2,M_grid,med_rel_mu(is,:),'-o', ...
        'LineWidth',1.5,'MarkerSize',4.5, ...
        'DisplayName',labels{is});
end
hold(ax2,'off');

set(ax2,'XScale','log','YScale','log');
xlabel(ax2,'Number of transitions, M');
ylabel(ax2,'Median relative error');
title(ax2,'(b) Recovery of \mu');
grid(ax2,'on');
box(ax2,'on');

%% (c) conditioning
ax3 = nexttile(tl,3);
hold(ax3,'on');

plot(ax3,c_grid,lambda_min,'-', ...
    'LineWidth',1.6, ...
    'DisplayName','\lambda_{min}(S^2)');

plot(ax3,scale_grid,lambda_sim,'o', ...
    'LineWidth',1.2, ...
    'MarkerSize',6, ...
    'DisplayName','Simulated scales');

hold(ax3,'off');

xlabel(ax3,'Mobility scale, c');
ylabel(ax3,'\lambda_{min}(S^2)');
title(ax3,'(c) Conditioning');
grid(ax3,'on');
box(ax3,'on');
legend(ax3,'Location','northeast','Box','off');

%% Common typography
axs = [ax1 ax2 ax3];
for k = 1:numel(axs)
    set(axs(k), ...
        'FontName','Times New Roman', ...
        'FontSize',9, ...
        'LineWidth',0.8, ...
        'TickDir','out');
end

%% Export
exportgraphics(fig,'experiment2_paper_figure_v2.pdf','ContentType','vector');
exportgraphics(fig,'experiment2_paper_figure_v2.png','Resolution',300);
savefig(fig,'experiment2_paper_figure_v2.fig');

fprintf('\nSaved:\n');
fprintf('  experiment2_paper_figure_v2.pdf\n');
fprintf('  experiment2_paper_figure_v2.png\n');
fprintf('  experiment2_paper_figure_v2.fig\n');

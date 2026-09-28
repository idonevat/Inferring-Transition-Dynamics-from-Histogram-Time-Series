%PLOT_T_SENSITIVITY Plot estimation quality versus observation length T.
%
% Main manuscript figure:
%   columns = source regime (equilibrium / nonequilibrium)
%   rows    = estimation error in pi, P, and mu
%
% Each panel compares the exact likelihood with the Gaussian approximation.
% The plotted statistic is the median estimation error across Monte Carlo
% replications.

clear; clc; close all;
S = load('experiment1_T_sensitivity_results.mat');

Tgrid = S.T_grid(:).';
nS = numel(S.scenario_name);
assert(nS==2,'Plotting script expects two source regimes.');

%% Main manuscript figure

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 24 22]);

tl = tiledlayout(fig,3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

lw = 1.5;
ms = 5.5;
metrics = {'err_pi','err_mu','err_P'};

ylabels = { ...
    'Median error in $\pi$', ...
    'Median relative error in $\mu$', ...
    'Median error in $P$'};


letters = 'abcdef';

for row = 1:3
    for s = 1:2

        panel = (row-1)*2+s;
        ax = nexttile(tl,panel);
        hold(ax,'on');

        fE = ['median_' metrics{row} '_exact'];
        fG = ['median_' metrics{row} '_gaussian'];

        yE = arrayfun(@(z) z.(fE),S.results.summary(s,:));
        yG = arrayfun(@(z) z.(fG),S.results.summary(s,:));

        plot(ax,Tgrid,yE,'-o', ...
            'LineWidth',lw, ...
            'MarkerSize',ms, ...
            'DisplayName','Exact likelihood');

        plot(ax,Tgrid,yG,'-s', ...
            'LineWidth',lw, ...
            'MarkerSize',ms, ...
            'DisplayName','Gaussian approximation');

        box(ax,'on');

        % Major grid only -- avoids clutter from log-scale minor grid lines
        grid(ax,'on');
        ax.XMinorGrid = 'off';
        ax.YMinorGrid = 'off';

        xlim(ax,[min(Tgrid) max(Tgrid)]);

        set(ax, ...
            'XScale','log', ...
            'XTick',Tgrid, ...
            'XTickLabel',string(Tgrid), ...
            'FontName','Times New Roman', ...
            'FontSize',13, ...
            'LineWidth',1, ...
            'TickDir','out');

        ylabel(ax,ylabels{row},'Interpreter','latex');

        % Only bottom row needs x-axis label
        if row == 3
            xlabel(ax,'Observation length, $T$', ...
                'Interpreter','latex');
        end

        % Clean panel labels
        text(ax,0.9,0.96,sprintf('(%c)',letters(panel)), ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'FontName','Times New Roman', ...
            'FontSize',14);

        % Legend constructed once
        if panel == 1
            lgd = legend(ax, ...
                'Location','northoutside', ...
                'Orientation','horizontal');
            lgd.Layout.Tile = 'north';
        end
    end
end

% Column headings
annotation(fig,'textbox',[0.2 0.935 0.22 0.035], ...
    'String','Equilibrium', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontName','Times New Roman', ...
    'FontSize',15);

annotation(fig,'textbox',[0.64 0.935 0.22 0.035], ...
    'String','Nonequilibrium', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontName','Times New Roman', ...
    'FontSize',15);

exportgraphics(fig, ...
    'experiment1_T_sensitivity.pdf', ...
    'ContentType','vector');

exportgraphics(fig, ...
    'experiment1_T_sensitivity.png', ...
    'Resolution',400);


%% Diagnostics table

fprintf('\nT-sensitivity summaries (median errors):\n');

for s = 1:nS

    fprintf('\n%s\n',S.scenario_name{s});

    fprintf([' T      pi-E     pi-G      P-E      P-G' ...
             '     mu-E     mu-G    badE%%   badG%%\n']);

    for iT = 1:numel(Tgrid)

        z = S.results.summary(s,iT);

        fprintf(['%3d  %8.4f %8.4f %8.4f %8.4f ' ...
                 '%8.4f %8.4f %7.2f %7.2f\n'], ...
            Tgrid(iT), ...
            z.median_err_pi_exact, ...
            z.median_err_pi_gaussian, ...
            z.median_err_P_exact, ...
            z.median_err_P_gaussian, ...
            z.median_err_mu_exact, ...
            z.median_err_mu_gaussian, ...
            100*z.frac_maxmu_gt10_exact, ...
            100*z.frac_maxmu_gt10_gaussian);

    end
end

fprintf('\nSaved T-sensitivity manuscript figure.\n');


% %PLOT_T_SENSITIVITY Plot estimation quality versus observation length T.
% %
% % Main figure: 3 rows (pi, P, mu) x 2 source regimes. Each panel compares
% % exact likelihood with the Gaussian approximation. The plotted statistic is
% % the median estimation error across Monte Carlo replications.
% 
% clear; clc; close all;
% S = load('experiment1_T_sensitivity_results.mat');
% 
% Tgrid = S.T_grid(:).';
% nS = numel(S.scenario_name);
% assert(nS==2,'Plotting script expects two source regimes.');
% 
% fig = figure('Color','w','Units','centimeters','Position',[2 2 24 24]);
% tl = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
% lw = 1.5; ms = 5.5;
% 
% metrics = {'err_pi','err_P','err_mu'};
% ylabels = {'Median error in $\pi$', ...
%            'Median error in $P$', ...
%            'Median relative error in $\mu$'};
% row_names = {'Equilibrium composition','Transition matrix','Interface mobility'};
% letters = 'abcdef';
% 
% for row=1:3
%     for s=1:2
%         ax = nexttile(tl,(row-1)*2+s);
%         hold(ax,'on');
%         fE = ['median_' metrics{row} '_exact'];
%         fG = ['median_' metrics{row} '_gaussian'];
%         yE = arrayfun(@(z) z.(fE),S.results.summary(s,:));
%         yG = arrayfun(@(z) z.(fG),S.results.summary(s,:));
% 
%         plot(ax,Tgrid,yE,'-o','LineWidth',lw,'MarkerSize',ms, ...
%             'DisplayName','Exact likelihood');
%         plot(ax,Tgrid,yG,'-s','LineWidth',lw,'MarkerSize',ms, ...
%             'DisplayName','Gaussian approximation');
%         grid(ax,'on'); box(ax,'on');
%         xlabel(ax,'Observation length, $T$','Interpreter','latex');
%         ylabel(ax,ylabels{row},'Interpreter','latex');
%         title(ax,sprintf('(%c) %s: %s',letters((row-1)*2+s), ...
%             row_names{row},S.scenario_name{s}),'FontWeight','normal');
%         xlim(ax,[min(Tgrid) max(Tgrid)]);
%         set(ax,'XScale','log','XTick',Tgrid,'XTickLabel',string(Tgrid), ...
%             'FontName','Times New Roman','FontSize',13,'LineWidth',1,'TickDir','out');
%         if row==1 && s==1
%             lgd = legend(ax,'Location','northoutside','Orientation','horizontal');
%             lgd.Layout.Tile = 'north';
%         end
%     end
% end
% 
% exportgraphics(fig,'experiment1_T_sensitivity_6panel.pdf','ContentType','vector');
% exportgraphics(fig,'experiment1_T_sensitivity_6panel.png','Resolution',400);
% 
% %% Compact pooled figure for manuscript use
% fig2 = figure('Color','w','Units','centimeters','Position',[2 2 23 9]);
% tl2 = tiledlayout(fig2,1,3,'TileSpacing','compact','Padding','compact');
% for row=1:3
%     ax = nexttile(tl2,row); hold(ax,'on');
%     % Pool both regimes at each T before taking the median.
%     yE = nan(size(Tgrid)); yG = nan(size(Tgrid));
%     for iT=1:numel(Tgrid)
%         xE = reshape(S.raw.([metrics{row} '_exact'])(:,iT,:),[],1);
%         xG = reshape(S.raw.([metrics{row} '_gaussian'])(:,iT,:),[],1);
%         yE(iT) = median(xE,'omitnan');
%         yG(iT) = median(xG,'omitnan');
%     end
%     plot(ax,Tgrid,yE,'-o','LineWidth',lw,'MarkerSize',ms,'DisplayName','Exact likelihood');
%     plot(ax,Tgrid,yG,'-s','LineWidth',lw,'MarkerSize',ms,'DisplayName','Gaussian approximation');
%     grid(ax,'on'); box(ax,'on');
%     xlabel(ax,'Observation length, $T$','Interpreter','latex');
%     ylabel(ax,ylabels{row},'Interpreter','latex');
%     title(ax,sprintf('(%c) %s','a'+row-1,row_names{row}),'FontWeight','normal');
%     set(ax,'XScale','log','XTick',Tgrid,'XTickLabel',string(Tgrid), ...
%         'FontName','Times New Roman','FontSize',13,'LineWidth',1,'TickDir','out');
%     if row==1
%         lgd = legend(ax,'Location','northoutside','Orientation','horizontal');
%         lgd.Layout.Tile = 'north';
%     end
% end
% exportgraphics(fig2,'experiment1_T_sensitivity_3panel.pdf','ContentType','vector');
% exportgraphics(fig2,'experiment1_T_sensitivity_3panel.png','Resolution',400);
% 
% %% Diagnostics table
% fprintf('\nT-sensitivity summaries (median errors):\n');
% for s=1:nS
%     fprintf('\n%s\n',S.scenario_name{s});
%     fprintf(' T      pi-E     pi-G      P-E      P-G     mu-E     mu-G    badE%%   badG%%\n');
%     for iT=1:numel(Tgrid)
%         z=S.results.summary(s,iT);
%         fprintf('%3d  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f %7.2f\n', ...
%             Tgrid(iT),z.median_err_pi_exact,z.median_err_pi_gaussian, ...
%             z.median_err_P_exact,z.median_err_P_gaussian, ...
%             z.median_err_mu_exact,z.median_err_mu_gaussian, ...
%             100*z.frac_maxmu_gt10_exact,100*z.frac_maxmu_gt10_gaussian);
%     end
% end
% 
% fprintf('\nSaved 6-panel diagnostic and 3-panel manuscript figures.\n');

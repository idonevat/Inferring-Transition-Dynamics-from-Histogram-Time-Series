%% plot_experiment1_QP_heatmaps.m
% Publication-style heatmaps for the true generator Q and one-step
% transition matrix P used in Experiment 1.
%
% Requires:
%   build_generator.m
%
% Produces:
%   experiment1_QP_heatmaps.pdf
%   experiment1_QP_heatmaps.png
%   experiment1_QP_heatmaps.fig

clear; clc;

%% True parameters used in Experiment 1
pi_true = [0.12 0.28 0.38 0.22];
mu_true = [0.35 0.25 0.30];
Delta = 1;

Q = build_generator(pi_true,mu_true);
P = expm(Delta*Q);

fprintf('True generator Q:\n');
disp(Q);

fprintf('One-step transition matrix P:\n');
disp(P);

%% Figure
fig = figure( ...
    'Units','inches', ...
    'Position',[1 1 7.2 3.2], ...
    'Color','w');

tl = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% Panel (a): Q
ax1 = nexttile(tl,1);
imagesc(ax1,Q);
axis(ax1,'image');
box(ax1,'on');

title(ax1,'(a) Generator Q');
xlabel(ax1,'Destination state');
ylabel(ax1,'Source state');

xticks(ax1,1:4);
yticks(ax1,1:4);

cb1 = colorbar(ax1);
cb1.Label.String = 'Transition rate';

% Add numerical values
for i = 1:size(Q,1)
    for j = 1:size(Q,2)
        text(ax1,j,i,sprintf('%.3f',Q(i,j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',13, ...
            'FontWeight','bold');
    end
end

%% Panel (b): P
ax2 = nexttile(tl,2);
imagesc(ax2,P);
axis(ax2,'image');
box(ax2,'on');

title(ax2,'(b) One-step transition matrix P');
xlabel(ax2,'Destination state');
ylabel(ax2,'Source state');

xticks(ax2,1:4);
yticks(ax2,1:4);

cb2 = colorbar(ax2);
cb2.Label.String = 'Transition probability';

% Add numerical values
for i = 1:size(P,1)
    for j = 1:size(P,2)
        text(ax2,j,i,sprintf('%.3f',P(i,j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',13, ...
            'FontWeight','bold');
    end
end

%% Typography
axs = [ax1 ax2];
for k = 1:numel(axs)
    set(axs(k), ...
        'FontName','Times New Roman', ...
        'FontSize',9, ...
        'LineWidth',0.8, ...
        'TickDir','out');
end

%% Export
exportgraphics(fig,'experiment1_QP_heatmaps.pdf','ContentType','vector');
exportgraphics(fig,'experiment1_QP_heatmaps.png','Resolution',300);
savefig(fig,'experiment1_QP_heatmaps.fig');

fprintf('\nSaved:\n');
fprintf('  experiment1_QP_heatmaps.pdf\n');
fprintf('  experiment1_QP_heatmaps.png\n');
fprintf('  experiment1_QP_heatmaps.fig\n');

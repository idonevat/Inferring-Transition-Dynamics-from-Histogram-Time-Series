%% plot_experiment4_final_figure.m
clear; clc;
repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
f1=fullfile(repoRoot,'results','experiment4','experiment4_stage1_population_results_R50.mat');
f2=fullfile(repoRoot,'results','experiment4','experiment4_stage2_temporal_results_R50.mat');
if ~isfile(f1), error('%s not found.',f1); end
if ~isfile(f2), error('%s not found.',f2); end
A=load(f1,'results'); pop=A.results; B=load(f2,'results'); tmp=B.results;
fig=figure('Color','w','Position',[60 80 1480 400]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax1=nexttile(tl,1); ik=numel(tmp.K_grid); iT=numel(tmp.T_grid);
mu_true=tmp.raw.true_mu{ik}; MU=tmp.raw.mu_hat{ik,iT}; MU=MU(all(isfinite(MU),2),:);
med=median(MU,1,'omitnan'); lo=prctile(MU,10,1); hi=prctile(MU,90,1); j=1:numel(mu_true);
fill(ax1,[j fliplr(j)],[lo fliplr(hi)],[0.88 0.88 0.88],'EdgeColor','none','FaceAlpha',0.55); hold(ax1,'on');
plot(ax1,j,mu_true,'k-','LineWidth',2.1); plot(ax1,j,med,'--','LineWidth',1.8); hold(ax1,'off'); box(ax1,'on');
xlabel(ax1,'State boundary j'); ylabel(ax1,'Transition intensity, \mu_j'); title(ax1,sprintf('(a) K=%d reconstruction',tmp.K_grid(ik)));
legend(ax1,{'10--90% Monte Carlo band','True profile','Median estimate'}, ...
    'Location','northeast','Box','off','FontSize',9.5); 
set(ax1,'FontSize',11,'LineWidth',0.8,'TickDir','out');

ax2=nexttile(tl,2); hold(ax2,'on');
for iK=1:numel(tmp.K_grid), plot(ax2,tmp.T_grid,tmp.med_err_mu(iK,:),'-o','LineWidth',1.8,'MarkerSize',6,'DisplayName',sprintf('K=%d',tmp.K_grid(iK))); end
refT=logspace(log10(min(tmp.T_grid)),log10(max(tmp.T_grid)),150); idx100=find(tmp.T_grid==100,1); anchor=median(tmp.med_err_mu(:,idx100),'omitnan'); ref=anchor*(refT/100).^(-0.5);
plot(ax2,refT,ref,'k--','LineWidth',1.35,'DisplayName','T^{-1/2} reference'); hold(ax2,'off'); box(ax2,'on'); grid(ax2,'on');
set(ax2,'XScale','log','YScale','log','FontSize',11,'LineWidth',0.8,'TickDir','out'); 
xlabel(ax2,'Number of observed transitions, T'); ylabel(ax2,'Median relative error in \mu'); 
title(ax2,'(b) Temporal resolution'); 
legend(ax2,'Location','southwest','Box','off','FontSize',9.5);

ax3=nexttile(tl,3); hold(ax3,'on');
for iO=1:numel(pop.occ_grid), plot(ax3,pop.K_grid,pop.med_time(:,iO),'-o','LineWidth',1.6,'MarkerSize',5,'DisplayName',sprintf('N/K=%d',pop.occ_grid(iO))); end
hold(ax3,'off'); box(ax3,'on'); grid(ax3,'on'); 
set(ax3,'XScale','log','YScale','log','FontSize',11,'LineWidth',0.8,'TickDir','out');
xlabel(ax3,'Number of states, K'); ylabel(ax3,'Median fitting time (s)'); 
title(ax3,'(c) Computational scaling'); 
legend(ax3,'Location','northwest','Box','off','FontSize',9.5);
base='experiment4_final_figure';
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
exportgraphics(fig,[base '.png'],'Resolution',400);
savefig(fig,[base '.fig']);

fprintf('Saved %s.pdf/png/fig\\n',base);
fprintf('\\nFINAL FIGURE SUMMARY\\n');
fprintf('K=%d,T=%d: median E_mu=%.4f, q90 E_mu=%.4f, median E_P=%.4f\\n',tmp.K_grid(end),tmp.T_grid(end),tmp.med_err_mu(end,end),tmp.q90_err_mu(end,end),tmp.med_err_P(end,end));
idx=max(1,numel(tmp.T_grid)-2):numel(tmp.T_grid); fprintf('Large-T slopes:\\n');
for iK=1:numel(tmp.K_grid), pp=polyfit(log(tmp.T_grid(idx)),log(tmp.med_err_mu(iK,idx)),1); fprintf(' K=%d: %.3f\\n',tmp.K_grid(iK),pp(1)); end
fprintf('Stage-1 median E_mu range: %.4f to %.4f\\n',min(pop.med_err_mu(:)),max(pop.med_err_mu(:)));
fprintf('Stage-1 median runtime at K=100: %s seconds\\n',mat2str(pop.med_time(end,:),4));

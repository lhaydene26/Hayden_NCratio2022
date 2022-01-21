% This script loads data generated with the compile_data.m and calculate_stats.m scripts
% and plots the nuclear spreading ratio for each genotype. This is the third and final script in the sequence of 
% quantifying nuclear spreading.


load('workspaces/compilation_ratios');
load('workspaces/compilation_ratios_binned');
load('workspaces/compilation_ratios_cc');


% for the specified cc (from s02_calculate_stats, run that first to generate dependencies), plot the statistics
% for the ratio
s = cell(size(compilation_ratios_cc,1),1);
fig = figure;
hold on;
conds = 1:5;
cmap = lines(numel(conds));
k = 1;
for i = conds
    x = k*ones(numel(compilation_ratios_cc{i,cc_num-3}(:,1)),1);
    x = x.*normrnd(1,1/60,numel(x),1);
    y = compilation_ratios_cc{i,cc_num-3}(:,1);
    
    scatter(x, y, 100,'filled','MarkerFaceColor',cmap(k,:),'LineWidth',10);
    err = errorbar(k,nanmean(y),nanstd(y),'color','k','LineWidth',4);
    err.Marker = '+';
    err.CapSize = 10;
    k = k + 1;
end
legend off
axis([0.5, numel(conds)+0.5, 0, 1]);
axis square
xticks(1:numel(conds))
ylabel('Length ratio (a.u.)','FontSize',36) % y-axis label
standardizePlot_bar_pvals(gcf,gca,'nuclearSpreadingRatios');
close(fig);


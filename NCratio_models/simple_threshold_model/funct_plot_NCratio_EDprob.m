function funct_plot_NCratio_EDprob(data_averaged, data_cc14_NCratio_avg, gs, data_fit_params)
    % This function plots probability of an embryo region dividing as a function of the 
    % N/C ratio, given a particular community radius
    
    x = find(gs==10);
    x(2) = find(gs==100);
    x(3) = find(gs==250);
    clear h
    clear xplot
    cmap = lines(numel(x));
    fig = figure;
    hold on
    xgrid = linspace(0, 1.5, 50);
    for i = 1:numel(x)
        xplot{i} = data_averaged(x(i),:,1)/data_cc14_NCratio_avg(x(i));
        plot(xplot{i}, data_averaged(x(i),:,2),'o','color',cmap(i,:),'LineWidth',2,'MarkerSize',15);
        h(i) = plot(xgrid, 1./(1+exp(-1.*(-xgrid-data_fit_params.pos(x(i)))./data_fit_params.steepness(x(i)))),'color',cmap(i,:),'LineWidth',2);
    end
    legend(h, '10 um','100', '250');
    legend boxoff
    xlabel('Relative N/C ratio')
    ylabel('Probability of extra division')
    standardizePlot(gcf,gca,strcat('figures/model_1/EDprob_vs_NCratio'));
    close(fig);

    
    
    
    
end

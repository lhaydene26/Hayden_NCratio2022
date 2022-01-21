function [] = standardizePlot(mygcf,mygca,printName)
    scaleFactor = 5;

    ax = mygca;
    ax.XLabel.FontSize = 7.5*scaleFactor;
    ax.YLabel.FontSize = 7.5*scaleFactor;
    ax.FontSize = 7.5*scaleFactor;
    set(mygca,'LineWidth',2,'FontName','Helvetica');
    set(mygcf,'units','centimeters');
    posgcf = mygcf.Position;
    posgca = mygca.Position;

    set(gcf,'units','centimeters','Position',[posgcf(1) posgcf(2) 5 5].*scaleFactor);
    set(gca,'Position',[0.19 0.16 0.7 0.7]);
    set(mygcf,'renderer','Painters')
    ax.Color = 'None';
    ax.XColor = 'black'; ax.YColor = 'black'; ax.ZColor = 'black';
    set(gcf,'PaperPositionMode','auto');
%     legend boxoff
    if(nargin==3)
        print('-depsc','-r300',printName,'-loose','-painters');
        saveas(gcf,strcat(printName,'.png'));
    end

end


function CNTPlotFast(fname,avgRate, steps, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad)
axesBound=12e-6;


xxx=numberBeams:-1:1;
not=ones(numberBeams,3).*NaN;

z(1,:)=1:numberBeams;
for bb=2:t
z(bb,:)=z(bb-1,:)+numberBeams;
end
zb=z(:); %List of nodes organized by CNT

plotmatrix=[nodeCoordinates(zb,1) nodeCoordinates(zb,2) nodeCoordinates(zb,3)]; %This needs to be in CNT down order

inserting=xxx*(t); %where to insert the NaN entries

slappy=insertrows(plotmatrix,not,[inserting]);

plot3(slappy(:,1)*10^6,slappy(:,2)*10^6,slappy(:,3)*10^6,'LineWidth',0.25,'Color','blue')


% plot(plotmatrix(:,:,1)*10^6, plotmatrix(:,:,2)*10^6,'LineWidth',0.5,'Color','blue')
%caxis([0, 10]);

    axis([0 5 0 5 0 10])
    figureHandle = gcf;
    set(gca,'FontSize',22)
    set(findall(gcf,'type','text'),'FontSize',22);
             set(gca,'FontSize',22)
             set(findall(gcf,'type','text'),'FontSize',22)
             %set(gca,'xtick',[])
             %set(gca,'XTick',[-2,0,2,4,6,8,10])
             xlabel('Substrate Position (\mum)');
             ylabel('Forest Height (\mum)');
             view(-184,5)
             %view(1)
             %view(2)
             
plotname=strcat(title,num2str(t));
saveas(gcf,fullfile(fname,plotname),'jpg');


close()
end
        


       
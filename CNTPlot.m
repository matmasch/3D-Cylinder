function CNTPlot(fname,avgRate, steps, h_span,numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad)
    
factor=1.2;    
%buffer=(22e-6-h_span)/2*1e6;
domain=0; %micron
%buffer=0;
 
%subplot(2,1,1)
for p=1:numberBeams; 
    pp=p:numberBeams:nodeCount;
    plot3(nodeCoordinates(pp,1)*1e6+domain,nodeCoordinates(pp,2)*1e6+domain, nodeCoordinates(pp,3)*1e6,'LineWidth',1);
    axis([-1 1 -1 1 0 2])
    figureHandle = gcf;
    set(gca,'FontSize',22)
    set(findall(gcf,'type','text'),'FontSize',22);
    %view([10,5,0])
    %xlabel('Substrate Position (\mum)');
    %ylabel('Forest Height (\mum)');
    hold on
    grid on
    %view(2)
end

       plotname=strcat(title,num2str(t));   
       saveas(gcf,fullfile(fname,plotname),'png');

            close()

% subplot(2,1,2)
%for p=1:numberBeams; 
%    pp=p:numberBeams:nodeCount;
%    plot3(nodeCoordinates(pp,1)*1e6+0.75,nodeCoordinates(pp,2)*1e6+0.75, nodeCoordinates(pp,3)*1e6,'LineWidth',0.5);
   
%    axis([0 3 0 3 0 6])
 %   view([10, 0, 0]);
%    figureHandle = gcf;
%    set(gca,'FontSize',22)
%    set(findall(gcf,'type','text'),'FontSize',22);
%    view([10,0,0])
%    xlabel('Substrate Position (\mum)');
%    ylabel('Forest Height (\mum)');
%    hold on
%end
%      plotname=strcat(title,'views',num2str(t));   
%       saveas(gcf,fullfile(fname,plotname),'png');
%      close()
% subplot(3,1,3)
% for p=1:numberBeams; 
%     pp=p:numberBeams:nodeCount;
%     plot3(nodeCoordinates(pp,1)*1e6+domain/2,nodeCoordinates(pp,2)*1e6+domain/2, nodeCoordinates(pp,3)*1e6,'LineWidth',0.5);
%    
%     axis([0 4 0 4 0 8])
%     view([0, 10, 0]);
%     figureHandle = gcf;
%     set(gca,'FontSize',22)
%     set(findall(gcf,'type','text'),'FontSize',22);
%     %view([10,5,0])
%     xlabel('Substrate Position (\mum)');
%     ylabel('Forest Height (\mum)');
%     hold on
% end

% if t==steps+totalCompress+1;
%    hold off
%    close()
%    j=steps+1:totalCompress+steps;
%    semilogy((j-steps)*0.05,-compressiveLoad(j))
%    figureHandle = gcf;
%    set(gca,'FontSize',16)
%    set(findall(gcf,'type','text'),'FontSize',16)
%    xlabel('Displacement (\mum)');
%    ylabel('Compressive Force (N)');       
% end
%         
             
%        plotname=strcat(title,num2str(t));
%              
%        saveas(gcf,fullfile(fname,plotname),'png');
       %saveas(gcf,fullfile(fname,plotname),'emf');
      %if ContinuousPlot==1 
%             close()
       %end
       
       
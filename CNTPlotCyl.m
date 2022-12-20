function CNTPlotCyl(fname,avgRate, steps, h_span,numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad)
    
factor=1.2;    
%buffer=(22e-6-h_span)/2*1e6;
domain=-0; %micron
%buffer=0;
 
% % subplot(2,2,1)
% % for p=1:numberBeams; 
% %     pp=p:numberBeams:nodeCount;
% %     plot3(nodeCoordinates(pp,1)*1e6+domain,nodeCoordinates(pp,2)*1e6+domain, nodeCoordinates(pp,3)*1e6,'LineWidth',0.1);
% %     axis([-10 10 -10 10 0 5])
% %     figureHandle = gcf;
% %     set(gca,'FontSize',22)
% %     set(findall(gcf,'type','text'),'FontSize',22);
% %     %view([10,5,0])
% %     %xlabel('Substrate Position (\mum)');
% %     %ylabel('Forest Height (\mum)');
% %     hold on
% %     grid on
% %     view(3)
% % end
% % subplot(2,2,2)
for p=1:numberBeams; 
    pp=p:numberBeams:nodeCount;
    plot3(nodeCoordinates(pp,1)*1e6+domain,nodeCoordinates(pp,2)*1e6+domain, nodeCoordinates(pp,3)*1e6,'LineWidth',0.1);
    axis([-5 5 -5 5 0 5])
    figureHandle = gcf;
    set(gca,'FontSize',22)
    set(findall(gcf,'type','text'),'FontSize',22);
    %view([10,5,0])
    %xlabel('Substrate Position (\mum)');
    %ylabel('Forest Height (\mum)');
    hold on
    grid on
    view(2)
end

       plotname=strcat(title,num2str(t));   
       saveas(gcf,fullfile(fname,plotname),'png');

            close()
  
       
       
       for i=1:numberBeams
       for j=i:numberBeams:nodeCount %this is the node numbers for
       nanotube(i,:)=j;
       end
       end
%        
%        CNTNodes(i,:)=nodeCoordinates(j,:);
% end
% end

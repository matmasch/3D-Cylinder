function CNTPlotBWfast(fname,avgRate, steps, h_span,numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,numberBeamsx,numberBeamsy)
    
count=0;  
%buffer=(22e-6-h_span)/2*1e6;
domain=-2; %micron
%buffer=0;
plotcomp=0
count=1;
for p=1:numberBeams;
    for pp=p:numberBeams:nodeCount;
        count=count+1;
        plotmatrix(count,:,:,:)=[nodeCoordinates(pp,1) nodeCoordinates(pp,2) nodeCoordinates(pp,3) ];
        plotcolor(count,:,:,:)=count/(1.2*numberBeamsx);
    end
       count=count+1;
       plotmatrix(count,:,:,:)=[NaN NaN NaN];
       plotcolor(count,:,:,:)=count/(1.2*numberBeamsx);
       end
plotcomp=1
enterningPlottingPart=1
    plot3(plotmatrix(:,1)*1e6,plotmatrix(:,2)*1e6, plotmatrix(:,3)*1e6,'LineWidth',1)%,'Color')%,plotcolor(:,:,:));
    axis([-4 4 -4 4 0 10])
    figureHandle = gcf;
    set(gca,'FontSize',22)
    set(findall(gcf,'type','text'),'FontSize',22);
    %view([10,5,0])
    %xlabel('Substrate Position (\mum)');
    %ylabel('Forest Height (\mum)');
    hold on
    grid on
    view(160,30)
       

% % %subplot(2,1,1)
% % for ff=1:numberBeamsy:numberBeams
% %     count=count+1;
% % for p=ff:ff+numberBeamsx; 
% %     pp=p:numberBeams:nodeCount;
% %     plot3(nodeCoordinates(pp,1)*1e6+domain,nodeCoordinates(pp,2)*1e6+domain, nodeCoordinates(pp,3)*1e6,'LineWidth',1,'Color',[0 0 0]+count/(1.2*numberBeamsx));
% %     axis([-4 4 -4 4 0 10])
% %     figureHandle = gcf;
% %     set(gca,'FontSize',22)
% %     set(findall(gcf,'type','text'),'FontSize',22);
% %     %view([10,5,0])
% %     %xlabel('Substrate Position (\mum)');
% %     %ylabel('Forest Height (\mum)');
% %     hold on
% %     grid on
% %     view(160,30)
% %     %view(2)
% % end
% % end

       plotname=strcat(title,num2str(t));   
       saveas(gcf,fullfile(fname,plotname),'jpg');

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
       
       
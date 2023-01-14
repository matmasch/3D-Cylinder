function CNTPlotBWPeriodic(fname,avgRate, steps, h_span,v_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,numberBeamsx,numberBeamsy,inactiveNodes)

allnodes=1:nodeCount;
activeNodes=setdiff(allnodes,inactiveNodes');

array = nodeCoordinates(activeNodes,1)/(1.05*h_span);
color = [0 0 0] + array;
scatter3(nodeCoordinates(activeNodes,1)*1e6,nodeCoordinates(activeNodes,2)*1e6, nodeCoordinates(activeNodes,3)*1e6,0.5,color, 'filled')
hold on         
axis off
            daspect([1 1 1])
            figureHandle = gcf;
            set(gca,'FontSize',22)
            hold on
            %view(96,15)
            view(96,5)
            
            
            % Add artificial points inbetween nodeCoordinates to make plots more
% full

l = find(activeNodes<nodeCount-numberBeams);
p=activeNodes(l);
fakeCoordinates(p,:)=(nodeCoordinates(p,:)+ nodeCoordinates((p+numberBeams),:))/2;
arrayf = fakeCoordinates(p,1)/(1.05*h_span);
colorf = [0 0 0] + arrayf;
scatter3(fakeCoordinates(p,1)*1e6,fakeCoordinates(p,2)*1e6, fakeCoordinates(p,3)*1e6,0.5,colorf, 'filled')

hold off
 
       plotname=strcat(fname,title,num2str(t));   
       print(gcf,plotname,'-dpng','-r300');
       %pause(5)
       %saveas(gcf,fullfile(fname,plotname),'jpg');
       %saveas(gcf,fullfile(fname,plotname),'fig');

            close()


       
       
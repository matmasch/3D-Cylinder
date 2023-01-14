function [closeNodes] =FindCloseNodesPar(nodeCoordinates,nodeCount,numberBeams)      
    gap=50e-9; % accidently 20 for IMECE study
    gapsquared=gap^2;
    i=[]; j=[]; closeNodes=[];
    

          parfor s=1:ceil( (nodeCount-numberBeams)*.6);
            for ss=s+1:ceil( (nodeCount-numberBeams)*.6);
                if abs(nodeCoordinates(s,1)-nodeCoordinates(ss,1))<gap;% & abs(nodeCoordinates(s,2)-nodeCoordinates(ss,2))<gap & abs(nodeCoordinates(s,3)-nodeCoordinates(ss,3))<gap;
                    sep=(nodeCoordinates(s,1)-nodeCoordinates(ss,1))^2 + (nodeCoordinates(s,2)-nodeCoordinates(ss,2))^2 + (nodeCoordinates(s,3)-nodeCoordinates(ss,3))^2;
                    if sep < gapsquared
                        closeNodes=[closeNodes; [s ss]];
                    end
 
                end
            end   
          end


          parfor s=ceil( (nodeCount-numberBeams)*.5): (nodeCount-numberBeams);
            for ss=s+1:(nodeCount-numberBeams);
                if abs(nodeCoordinates(s,1)-nodeCoordinates(ss,1))<gap;% & abs(nodeCoordinates(s,2)-nodeCoordinates(ss,2))<gap & abs(nodeCoordinates(s,3)-nodeCoordinates(ss,3))<gap;
                    sep=(nodeCoordinates(s,1)-nodeCoordinates(ss,1))^2 + (nodeCoordinates(s,2)-nodeCoordinates(ss,2))^2 + (nodeCoordinates(s,3)-nodeCoordinates(ss,3))^2;
                    if sep < gapsquared
                        closeNodes=[closeNodes; [s ss]];
                    end
 
                end
            end   
          end


    
    if size(closeNodes,1)==0;
        closeNodes=0;
    end
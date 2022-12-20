if size(closeNodes,1)==0; 
         closeNodes=0; %%[nodeCount,0] closeNodes==[nodeCount,0]
     else
         sizeClose=size(closeNodes,1);
         xa=zeros(sizeClose,1); ya=zeros(sizeClose,2); za=zeros(sizeClose,3); LL=zeros(sizeClose,1);
         CC=zeros(sizeClose,1); SS=zeros(sizeClose,1);
         xa=nodeCoordinates(closeNodes(:,2),1)-nodeCoordinates(closeNodes(:,1),1);
         ya=nodeCoordinates(closeNodes(:,2),2)-nodeCoordinates(closeNodes(:,1),2);
         za=nodeCoordinates(closeNodes(:,2),3)-nodeCoordinates(closeNodes(:,1),3);
             if PeriodicBoundary>0;
                for hh=1:sizeClose;
                   if abs(xa)>h_span/2;
                        if nodeCoordinates(closeNodes(hh,1),1)>nodeCoordinates(closeNodes(hh,2),1) 
                            xa=span+nodeCoordinates(closeNodes(hh,2),1)-nodeCoordinates(closeNodes(h,1),1);
                            else
                                xa=span-nodeCoordinates(closeNodes(hh,2),1)+nodeCoordinates(closeNodes(hh,1),1);
                            end
                    end
                end
             end
         
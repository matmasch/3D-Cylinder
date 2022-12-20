function [gap, sep,closeNodes] =FindCloseNodesSparseCyl2(nodeCoordinates,nodeCount,numberBeamsx,numberBeamsy,t,rmax)      
 
    sep=sparse(1,1,0,nodeCount,nodeCount);
    gap=50e-9;
    
    if rmax < 15e-9;
        gap=50e-9;
    end
    %gap=60e-9;%%
    counter=1;
    jj(1)=1;  ii(1)=2; dist(1)=1;
    numberBeamsRad=numberBeamsx;     numberBeams=numberBeamsx*numberBeamsy;
    CNT=zeros(numberBeamsx/2, numberBeamsy);
    
    
        i=1:numberBeamsRad/2;    
        for j=1:numberBeamsy;
        CNT(:,j)=i+numberBeamsRad*(j-1);
        end
        CNT=CNT(:);
        list=zeros(numberBeamsx/2*numberBeamsy,t);

        for i=1:size(CNT,1);
            list(i,:)=CNT(i):numberBeams:numberBeams*t;
        end
        list=list(:);
            length=size(list,1); %number of elements in list
    
    

for s=1:length;
            for ss=s+1:length;
                if abs(nodeCoordinates(list(s),1)-nodeCoordinates(list(ss),1))<gap && abs(nodeCoordinates(list(s),2)-nodeCoordinates(list(ss),2))<gap ;
                    distance=sqrt(   (nodeCoordinates(list(s),1)-nodeCoordinates(list(ss),1))^2 + (nodeCoordinates(list(s),2)-nodeCoordinates(list(ss),2))^2+ (nodeCoordinates(list(s),3)-nodeCoordinates(list(ss),3))^2  );
                    if distance<gap
                    jj(counter)=list(ss);
                    ii(counter)=list(s);
                    dist(counter)=distance;
                    counter=counter+1;
                   end
                end
            end
end






CNT=zeros(numberBeamsx/2, numberBeamsy);
        i=numberBeamsRad/2+1:numberBeamsRad;    
        for j=1:numberBeamsy;
        CNT(:,j)=i+numberBeamsRad*(j-1);
        end
        CNT=CNT(:);
        list=zeros(numberBeamsx/2*numberBeamsy,t);

        for i=1:size(CNT,1);
            list(i,:)=CNT(i):numberBeams:numberBeams*t;
        end
        list=list(:);
            length=size(list,1); %number of elements in list

for s=1:length;
            for ss=s+1:length;
                if abs(nodeCoordinates(list(s),1)-nodeCoordinates(list(ss),1))<gap && abs(nodeCoordinates(list(s),2)-nodeCoordinates(list(ss),2))<gap ;
                    distance=sqrt(   (nodeCoordinates(list(s),1)-nodeCoordinates(list(ss),1))^2 + (nodeCoordinates(list(s),2)-nodeCoordinates(list(ss),2))^2+ (nodeCoordinates(list(s),3)-nodeCoordinates(list(ss),3))^2  );
                    if distance<gap
                    jj(counter)=list(ss);
                    ii(counter)=list(s);
                    dist(counter)=distance;
                    counter=counter+1;
                   end
                end
            end
end

sep=sparse(jj(:), ii(:), dist(:), nodeCount, nodeCount);

    closeNodes=[ii(:),jj(:)];
    closeNodes=unique(closeNodes,'rows');
    if counter==1;
        closeNodes=0;
    end

    
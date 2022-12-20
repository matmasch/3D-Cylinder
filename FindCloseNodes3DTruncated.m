function [gap, sep,closeNodes] =FindCloseNodes3DTruncated(nodeCoordinates,nodeCount,numberBeamsx,numberBeamsy,t,steps,oldCloseNodes,truncatedTime)      
 
    sep=sparse(1,1,0,nodeCount,nodeCount);
    gap=50e-9;
    numberBeams=numberBeamsx*numberBeamsy;
    
    %truncatedTime=50;
    truncation=truncatedTime*numberBeams;
    
%    if rmax < 15e-9;
        gap=50e-9;
   % end
    %gap=60e-9;%%
    counter=1; ii=0; jj=0;
    dist=0; jj(1)=1;  ii(1)=2; dist(1)=1;
    
    
    
 if t<truncatedTime+1   
 list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
 CNT=zeros(numberBeamsx/2, numberBeamsy/2);

 i=1:numberBeamsx/2;    
for j=1:numberBeamsy/2;
CNT(:,j)=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,t-1);
for i=1:size(CNT,1);
    list(i,:)=CNT(i):numberBeams:numberBeams*(t-1);
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






list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
CNT=zeros(numberBeamsx/2, numberBeamsy/2);

i=1:numberBeamsx/2; 
for j=numberBeamsy/2+1:numberBeamsy;
    CNT(:,(j-numberBeamsy/2))=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,t-1);
for i=1:size(CNT,1);
    list(i,:)=CNT(i):numberBeams:numberBeams*(t-1);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
 CNT=zeros(numberBeamsx/2, numberBeamsy/2);
i=numberBeamsx/2+1:numberBeamsx;    

for j=numberBeamsy/2+1:numberBeamsy;
    CNT(:,(j-numberBeamsy/2))=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,t-1);
for i=1:size(CNT,1);
    list(i,:)=CNT(i):numberBeams:numberBeams*(t-1);
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


list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
CNT=zeros(numberBeamsx/2, numberBeamsy/2);

 i=numberBeamsx/2+1:numberBeamsx;    
for j=1:numberBeamsy/2;
    CNT(:,j)=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,t-1);
for i=1:size(CNT,1);
    list(i,:)=CNT(i):numberBeams:numberBeams*(t-1);
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




   % contact=find(sep < gap & sep>0);%%Returns the node numbers of CNTs within a given distance

    %[II,JJ]=ind2sub(size(sep),contact);%%converts from index to column,row notation
    closeNodes=[ii(:),jj(:)];
    closeNodes=unique(closeNodes,'rows');
    if counter==1;
        closeNodes=0;
    end

 end  
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nodeCount>truncation
    sep=sparse(1,2,1,truncation,truncation);
    nodeUpperLimit=truncation;
    start=nodeCount-nodeUpperLimit;
    
 CNT=zeros(numberBeamsx/2, numberBeamsy/2);

 i=1:numberBeamsx/2;    
for j=1:numberBeamsy/2;
CNT(:,j)=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,truncatedTime-2);
for i=1:size(CNT,1);
    list(i,:)=(CNT(i)+start):numberBeams:numberBeams*(t-1);
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






list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
CNT=zeros(numberBeamsx/2, numberBeamsy/2);

i=1:numberBeamsx/2; 
for j=numberBeamsy/2+1:numberBeamsy;
    CNT(:,(j-numberBeamsy/2))=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,truncatedTime-2);
for i=1:size(CNT,1);
    list(i,:)=CNT(i)+start:numberBeams:numberBeams*(t-1);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
 CNT=zeros(numberBeamsx/2, numberBeamsy/2);
i=numberBeamsx/2+1:numberBeamsx;    

for j=numberBeamsy/2+1:numberBeamsy;
    CNT(:,(j-numberBeamsy/2))=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,truncatedTime-2);
for i=1:size(CNT,1);
    list(i,:)=CNT(i)+start:numberBeams:numberBeams*(t-1);
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


list=zeros(numberBeamsx/2, numberBeamsy/2); %Splits domain into quadrant
CNT=zeros(numberBeamsx/2, numberBeamsy/2);

 i=numberBeamsx/2+1:numberBeamsx;    
for j=1:numberBeamsy/2;
    CNT(:,j)=i+numberBeamsx*(j-1);
end
CNT=CNT(:);

list=zeros(numberBeamsx/2*numberBeamsy/2,truncatedTime-2);
for i=1:size(CNT,1);
    list(i,:)=CNT(i)+start:numberBeams:numberBeams*(t-1);
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

    
        if size(closeNodes,2)>1
        closeNodes=[closeNodes; oldCloseNodes];
    end
    if counter==1;
        closeNodes=0;
    end
    closeNodes=unique(closeNodes,'rows');

end
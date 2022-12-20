function [gap, sep,closeNodes] =FindCloseNodesSparse1(nodeCoordinates,nodeCount,numberBeams,t,rmax)      
 
    sep=sparse(1,1,0,nodeCount,nodeCount);
    gap=50e-9;
    
    if rmax < 15e-9;
        gap=50e-9;
    end
    %gap=60e-9;%%
    counter=1;
    jj(1)=1;
    ii(1)=2;
    dist(1)=1;

if t<30;
        for s=1:(nodeCount-numberBeams);
            for ss=s+1:(nodeCount-numberBeams);
                if abs(nodeCoordinates(s,1)-nodeCoordinates(ss,1))<gap && abs(nodeCoordinates(s,2)-nodeCoordinates(ss,2))<gap ;
                    dist(counter)=sqrt(   (nodeCoordinates(s,1)-nodeCoordinates(ss,1))^2 + (nodeCoordinates(s,2)-nodeCoordinates(ss,2))^2+ (nodeCoordinates(s,3)-nodeCoordinates(ss,3))^2  );
                    jj(counter)=ss;
                    ii(counter)=s;
                    if dist(counter)<gap
                        counter=counter+1;
                    end
 
                end   
            end
        end
    sep=sparse(jj(:), ii(:), dist(:), nodeCount, nodeCount);
end

if t>29 
       for s=1:(nodeCount-numberBeams)/2;
        for ss = s+1:(nodeCount-numberBeams)/2;
           if ((nodeCoordinates(s,1)-nodeCoordinates(ss,1)) <= gap &&  (nodeCoordinates(s,2)-nodeCoordinates(ss,2)) <= gap);
                    jj(counter)=ss;
                    ii(counter)=s;
                    dist(counter)=sqrt(   (nodeCoordinates(s,1)-nodeCoordinates(ss,1))^2 + (nodeCoordinates(s,2)-nodeCoordinates(ss,2))^2+ (nodeCoordinates(s,3)-nodeCoordinates(ss,3))^2  );
                    if dist(counter)<gap
                        counter=counter+1;
                    end
           end
        end
       end

  
       for s=floor((nodeCount-numberBeams)*1/4):floor((nodeCount-numberBeams)*3/4);
        for ss=max(s+1,(nodeCount-numberBeams)/2) :floor((nodeCount-numberBeams)*3/4);
             if ((nodeCoordinates(s,1)-nodeCoordinates(ss,1)) <= gap &&  (nodeCoordinates(s,2)-nodeCoordinates(ss,2)) <= gap);
                                jj(counter)=ss;
                    ii(counter)=s;
                    dist(counter)=sqrt(   (nodeCoordinates(s,1)-nodeCoordinates(ss,1))^2 + (nodeCoordinates(s,2)-nodeCoordinates(ss,2))^2+ (nodeCoordinates(s,3)-nodeCoordinates(ss,3))^2  );
                    if dist(counter)<gap
                        counter=counter+1;
                    end
             end
       end
       end
    
           for  s=(nodeCount-numberBeams)/2:(nodeCount-numberBeams);
                for ss=max(s+1,floor((nodeCount-numberBeams)*.75)):(nodeCount-numberBeams);
                    if ((nodeCoordinates(s,1)-nodeCoordinates(ss,1)) <= gap &&  (nodeCoordinates(s,2)-nodeCoordinates(ss,2)) <= gap);
                                jj(counter)=ss;
                    ii(counter)=s;
                    dist(counter)=sqrt(   (nodeCoordinates(s,1)-nodeCoordinates(ss,1))^2 + (nodeCoordinates(s,2)-nodeCoordinates(ss,2))^2+ (nodeCoordinates(s,3)-nodeCoordinates(ss,3))^2  );
                    if dist(counter)<gap
                        counter=counter+1;
                    end
                end
           end
           end
    sep=sparse(jj(:), ii(:), dist(:), nodeCount, nodeCount);
end

     
    
   % contact=find(sep < gap & sep>0);%%Returns the node numbers of CNTs within a given distance

    %[II,JJ]=ind2sub(size(sep),contact);%%converts from index to column,row notation
    closeNodes=[ii(:),jj(:)];
    closeNodes=unique(closeNodes,'rows');
    if counter==1;
        closeNodes=0;
    end

    
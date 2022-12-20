function [elementNodes,ang,rate,nodeCoordinates,nucleationSite,growthNode,ro,phi]= nucleateRadial2(numberBeamsx,numberBeamsy,diam,length,avgRate,rout,rate_stdev,ang_stdev);

radius=diam/2*10^-6;
numberBeamsRad=numberBeamsx;
numberBeams=numberBeamsx*numberBeamsy;

for o=1:numberBeamsRad %  Indexes at a constant angle first
    for count=1:numberBeamsy % Then indexes along the length
    % uniformly spaced nucleation points along 1/4 perimeter and allows buffer of 3 CNTs at edges
        index=(count-1)*numberBeamsx+o;
        origin(index,1)=radius*cos(-pi/2+2*pi/(numberBeamsRad-1)*(o-1)); %x-Direction
        origin(index,2)=radius*sin(-pi/2+2*pi/(numberBeamsRad-1)*(o-1)); %y-Direction
        origin(index,3)=length*10^-6*(count-1)/numberBeamsy; %Length is in the z-direction
        phi(index)= pi/2+ .2*rand;
    end
end

   for ii=1:numberBeams;
     mu_ang = [0];
     Sigma_ang = [(ang_stdev*pi/180)^2;]; 
     R_ang = chol(Sigma_ang);
     rand_ang(ii) = repmat(mu_ang,1) + randn(1)*R_ang;
     %phi(ii)=pi/2+ (rand-0.5)*rand_ang(ii); %Used to be ang(ii)=
     ro(ii)=rout;%%30e-9+8e-9*rand()
     
     mu_rate=avgRate;
     Sigma_rate = [(rate_stdev/100*avgRate)^2;];
     R_rate = chol(Sigma_rate);
     rate(ii) = repmat(mu_rate,1) + randn(1)*R_rate;  
   end


for a=1:numberBeamsRad/2;
    for counter=1:numberBeamsy  
    cc=(counter-1)*numberBeamsRad+a;
    omega=atan(origin(cc,2)/origin(cc,1)); %%CHANGED THIS - Was 2/1
    ang(cc)=omega+rand_ang(cc);%% this is the global growth angle
end
end

for a=numberBeamsRad/2+1:numberBeamsRad;
for counter=1:numberBeamsy
    cc=(counter-1)*numberBeamsRad+a;
    omega=atan(origin(cc,2)/origin(cc,1))+pi; %% CHANGED THIS - WAS 2/1
    ang(cc)=omega+rand_ang(cc);%% this is the global growth angle
end
end 
 


   
   
nodeCount=0;
element=0;

numberBeams=numberBeamsx*numberBeamsy;
for num = numberBeams+1:2*numberBeams; %% setting up CNT bases
    nodeCount=nodeCount+1;
    xxx(num)=origin(num-numberBeams,1);
    yyy(num)=origin(num-numberBeams,2); 
    zzz(num)=origin(num-numberBeams,3);
    nodeCoordinates(num,1)=xxx(num);
    nucleationSite(num-numberBeams,1)=xxx(num);
    nodeCoordinates(num,2)=yyy(num);
    nucleationSite(num-numberBeams,2)=yyy(num);
    nodeCoordinates(num,3)=zzz(num);
    nucleationSite(num-numberBeams,3)=zzz(num);
end



 %% making the free tips of the CNTs the initial nodes
    
     
         for ii=1:numberBeamsy
             for num=(ii-1)*numberBeamsRad+1:(ii-1)*numberBeamsRad+numberBeamsRad/2;%%Setting position of CNT free ends
        element=element+1; 
        nodeCount=nodeCount+1;
        %xxx(num)=xxx(num+numberBeams)+cos(ang(num))*rate(num);
        xxx(num)=xxx(num+numberBeams)+cos(ang(num))*sin(phi(num))*rate(num);
        growthNode(num,1)=xxx(num);
        %yyy(num)=yyy(num+numberBeams)+(sin(ang(num))*rate(num));
        yyy(num)=yyy(num+numberBeams)+sin(ang(num))*sin(phi(num))*rate(num);
        growthNode(num,2)=yyy(num);
        %zzz(num)=zzz(num+numberBeams)+(sin(ang(num))*rate(num));
        zzz(num)=zzz(num+numberBeams)+cos(phi(num))*rate(num);
        growthNode(num,3)=zzz(num);
        
        
% % % %     nodeCoordinates(num,1)=nodeCoordinates(num+numberBeams,1)+sin(ang(num))*cos(phi(num))*rate(num) %changed '*' to '+' before sin
% % % %     nodeCoordinates(num,2)=nodeCoordinates(num+numberBeams,2)+sin(ang(num))*sin(phi(num))*rate(num) %changed '*' to '+' before sin
% % % %     nodeCoordinates(num,3)=nodeCoordinates(num+numberBeams,3)+cos(ang(num))*rate(num) %changed '*' to '+' before cos
        
        
        elementNodes(element,1)=element;
        elementNodes(element,2)=nodeCount;
        
        nodeCoordinates(num,1)=xxx(num);
        nodeCoordinates(num,2)=yyy(num);
        nodeCoordinates(num,3)=zzz(num);
     end  
    end
       
   for ii=1:numberBeamsy
     for num=(ii-1)*numberBeamsRad+numberBeamsRad/2+1:(ii-1)*numberBeamsRad+numberBeamsRad;%%Setting position of CNT free ends
          
        element=element+1; 
        nodeCount=nodeCount+1;
        %xxx(num)=xxx(num+numberBeams)+cos(ang(num))*rate(num);
        xxx(num)=xxx(num+numberBeams)+cos(ang(num))*sin(phi(num))*rate(num);
        growthNode(num,1)=xxx(num);
        %yyy(num)=yyy(num+numberBeams)+(sin(ang(num))*rate(num));
        yyy(num)=yyy(num+numberBeams)+sin(ang(num))*sin(phi(num))*rate(num);
        growthNode(num,2)=yyy(num);
        %zzz(num)=zzz(num+numberBeams)+(sin(ang(num))*rate(num));
        zzz(num)=zzz(num+numberBeams)+cos(phi(num))*rate(num);
        growthNode(num,3)=zzz(num);
        
        elementNodes(element,1)=element;
        elementNodes(element,2)=nodeCount;
        
        nodeCoordinates(num,1)=xxx(num);
        nodeCoordinates(num,2)=yyy(num);
        nodeCoordinates(num,3)=zzz(num);
     end  
    end
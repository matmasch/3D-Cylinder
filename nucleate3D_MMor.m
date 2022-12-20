%function [elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,ro]=nucleate10umCompress(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rout,rate_stdev,ang_stdev)

h_span=5;
v_span=5;
numberBeamsx=10;
numberBeamsy=10;

ang_stdev=10;
rate_stdev=5;
avgRate=50e-9;
rout=10e-9;


for o=1:numberBeamsx
    o
    origin1(o)=h_span/(numberBeamsx)*(o-1/2)
end

for p=1:numberBeamsy
    origin2(p)=v_span/(numberBeamsy)*(p-1/2)
end
    [x,y]=meshgrid(origin1,origin2)
    numberBeams=numberBeamsx*numberBeamsy
    
    
    
for ii=1:numberBeams;
    mu_ang = 0;
    Sigma_ang = (ang_stdev*3.1415/180)^2; %Converted to radians
    R_ang = chol(Sigma_ang);
    ang(ii) = repmat(mu_ang,1) + randn(1)*R_ang;
    
    ro(ii)=rout;
    phi(ii)=rand*pi;
    
    mu_rate=avgRate;
    Sigma_rate = [(rate_stdev/100*avgRate)^2;];%changed
    R_rate = chol(Sigma_rate);
    rate(ii) = repmat(mu_rate,1) + randn(1)*R_rate;  
end



nodeCount=0;
element=0;

%% Setting up CNT bases 
for num = numberBeams+1:2*numberBeams
    num
    nodeCount=nodeCount+1
    nodeCoordinates(num,1)=x(num-numberBeams)
    nodeCoordinates(num,2)=y(num-numberBeams) 
    nodeCoordinates(num,3)=0
%%Making the free tips of the CNTs the initial nodes

%% Matt Note:  I am not sure what these lines of code are doing
    %nucleationSite(num-numberBeams,1)=nodeCoordinates(num,1);%%xxx(num)
    %nucleationSite(num-numberBeams,2)=nodeCoordinates(num,2);%%yyy(num)
    %nucleationSite(num-numberBeams,3)=nodeCoordinates(num,3);%%zzz(num)
end


for num=1:numberBeams %%Setting position of CNT free ends
    num
    element=element+1 
    nodeCount=nodeCount+1
    nodeCoordinates(num,1)=nodeCoordinates(num+numberBeams,1)+sin(ang(num))*cos(phi(num))*rate(num); %changed '*' to '+' before sin
    nodeCoordinates(num,2)=nodeCoordinates(num+numberBeams,2)+sin(ang(num))*sin(phi(num))*rate(num); %changed '*' to '+' before sin
    nodeCoordinates(num,3)=nodeCoordinates(num+numberBeams,3)+cos(ang(num))*rate(num) %changed '*' to '+' before cos
    growthNode(num,1)=nodeCoordinates(num,1)
    growthNode(num,2)=nodeCoordinates(num,2)
    elementNodes(element,1)=element
    elementNodes(element,2)=nodeCount
end   

         %%Adding Plot Function%%
     
             for p=1:numberBeams
                 p
                 pp=p:numberBeams:nodeCount
                 %plot(nodeCoordinates(pp,1)*1e6+buffer*1e6,nodeCoordinates(pp,2)*1e6,'Marker','o','MarkerSize',4,'LineWidth',1);
                 plot3(nodeCoordinates(pp,1)*1e6,nodeCoordinates(pp,2)*1e6, nodeCoordinates(pp,3)*1e6,'LineWidth',1)
                 %plot(nodeCoordinates(pp,1),nodeCoordinates(pp,2),'LineWidth',1);
                 %scatter(nodeCoordinates(closeNodes(:,1),1),nodeCoordinates(closeNodes(:,1),2))
                 %axis([0 axesBound*1e6 -2 10]);
                 figureHandle = gcf
                 set(gca,'FontSize',22)
                 set(findall(gcf,'type','text'),'FontSize',22)
                 %set(gca,'xtick',[])
                 set(gca,'XTick',[-2,0,2,4,6,8,10])
                 xlabel('Substrate Position (\mum)')
                 ylabel('Forest Height (\mum)')
                 hold on
             end


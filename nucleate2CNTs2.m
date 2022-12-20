function [elementNodes,ang,rate,nodeCoordinates,nucleationSite,growthNode,ro,phi]= nucleate2CNTs2(rout)

numberBeams=2;

x(1)=0; x(2)=250e-9;
y(1)=0; y(2)=0;
z(1)=0; z(2)=0;

ang(1)=0.01;     ang(2)=0.2474; %Angle with z-axis
phi(1)=0;        phi(2)=pi; %X-Y plane
rate(1)=50e-9;   rate(2)=55e-9;
ro(1)=rout;      ro(2)=rout;

nodeCount=0; element=0;

for num = numberBeams+1:2*numberBeams %% Setting up CNT bases
    nodeCount=nodeCount+1;
    nodeCoordinates(num,1)=x(num-numberBeams);
    nodeCoordinates(num,2)=y(num-numberBeams); 
    nodeCoordinates(num,3)=0;
    nucleationSite(num-numberBeams,1)=nodeCoordinates(num,1);%%xxx(num)
    nucleationSite(num-numberBeams,2)=nodeCoordinates(num,2);%%yyy(num)
    nucleationSite(num-numberBeams,3)=nodeCoordinates(num,3);%%zzz(num)
end

%%Making the free tips of the CNTs the initial nodes
for num=1:numberBeams;%%Setting position of CNT free ends
        element=element+1; 
        nodeCount=nodeCount+1;
        nodeCoordinates(num,1)=nodeCoordinates(num+numberBeams,1)+sin(ang(num))*cos(phi(num))*rate(num);%changed '*' to '+' before sin
        nodeCoordinates(num,2)=nodeCoordinates(num+numberBeams,2)+sin(ang(num))*sin(phi(num))*rate(num);%changed '*' to '+' before sin
        nodeCoordinates(num,3)=nodeCoordinates(num+numberBeams,3)+cos(ang(num))*rate(num);%changed '*' to '+' before cos
        growthNode(num,1)=nodeCoordinates(num,1);
        growthNode(num,2)=nodeCoordinates(num,2);
        growthNode(num,3)=nodeCoordinates(num,3);
        elementNodes(element,1)=element;
        elementNodes(element,2)=nodeCount;
end   
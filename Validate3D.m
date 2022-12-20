clear all
E=1E12; ro=5E-9; ri=3E-9; G=5e11;

nodeCoordinates=[0.25E-6 1E-6 0;0 0 0]; % Node Coordinates
elementNodes=[1 2]; %Connectivity matrix
element=1;

A=pi*(ro.^2-ri.^2)'; Ic=pi/4*(ro.^4-ri.^4)'; GDof=12;

    xa=nodeCoordinates(elementNodes(:,1),1)-nodeCoordinates(elementNodes(:,2),1);
    ya=nodeCoordinates(elementNodes(:,1),2)-nodeCoordinates(elementNodes(:,2),2);
    za=nodeCoordinates(elementNodes(:,1),3)-nodeCoordinates(elementNodes(:,2),3);

    LL=(xa.*xa+ ya.*ya + za.*za).^0.5;
    l=xa./LL;     m=ya./LL;     n=za./LL; 

    D=sqrt(l.*l+m.*m);
    Lvx=-m./D;     Lvy=l./D;     Lvz=0./D;
    %Lvector2=[Lvx Lvy Lvx];
    
    Lv3x=-(l.*n./D);
    Lv3y=-(m.*n./D);
    
CXx = (xa)./LL; CYx = (ya)./LL; CZx = (za)./LL;
D = sqrt(CXx.*CXx + CYx.*CYx);
CXy = -CYx./D;
CYy = CXx./D;
CZy = 0;
CXz = -CXx.*CZx./D;
CYz = -CYx.*CZx./D;
CZz = D;

K= SpaceBeamElemStiffness_mod2(E,A,GDof,LL,Ic,G,elementNodes,element,CXx,CYx,CZx,CXy,CYy,CZy,CXz,CYz,CZz);

force=zeros(12,1);
force(2)=-5E-9

alldof=[7:12];
prescribedDof=[alldof]';
u=solution2(GDof,prescribedDof,K,force)

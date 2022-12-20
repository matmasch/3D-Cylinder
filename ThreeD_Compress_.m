E=1E12; rout=12.5E-9; ri=7.5E-9; A=pi*(rout^2-ri^2); I=pi/4*(rout^4-ri^4); EA=E*A; EI=E*I;
G=5e11;

load('C:\Users\MattMaschmann\Documents\MATLAB\ImageFolder\Official\Stats\150_10um1\Grow150CNT_10um_test1_.mat')
fname = 'C:\Users\MattMaschmann\Documents\MATLAB\ImageFolder\Temperature\ThemoMechCompress';
title=input('Name of Run ','s')

totalCompress=130;
compressiveLoad=0; %% Initialization Only
%%radius=2.5e-6;%% used for radial systems
conductivity=1000;% CNT conductivity
ContinuousPlot=1;%% 0=plotting off.  1=plotting on
PeriodicBoundary=0; %% 0=off, 1=on
%beamType=0; %%0=Euler Beam ; 1=Timoshenko Beam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              COMPRESSION OF CNT FOREST GROWN ABOVE                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GDof=6*size(nodeCoordinates,1);
UU=zeros(GDof,1); forceCompress=zeros(GDof,1); U=zeros(GDof,1); compressiveLoad=zeros(totalCompress,1);

for e=element-numberBeams+1:element;
        L(e)=(rate(e-element+numberBeams));
        A(e)=(A(e-element+numberBeams));
        I(e)=I(e-element+numberBeams)';
end

topNode=max(nodeCoordinates(:,3)); %Top-most node in forest


%% COMPRESSION LOOP %%%%%
for compress=1:totalCompress;
    U=zeros(GDof,1); forceCompress=zeros(GDof,1); UU=sparse(GDof,1,0,GDof,1); activeForce=zeros(GDof,1);
    K=sparse(5,5,0,GDof,GDof);  AA=sparse(5,5,0,GDof,GDof); 
      
    compress
    t=compress+steps;
    [UU,compressedNodes] = nodeCompression3D(topNode,nodeCoordinates,compress,GDof);
    [gap, sep,closeNodes] =FindCloseNodesSparse2(nodeCoordinates,nodeCount,numberBeamsx,numberBeamsy,t,rmax)  ;
    %[closeNodes]=FindCloseNodesSparse2(nodeCoordinates,nodeCount,numberBeams,t,rmax,leftTrim,trimmed);
 
% % % % % % %     xa=nodeCoordinates(elementNodes(:,2),1)-nodeCoordinates(elementNodes(:,1),1);
% % % % % % %     ya=nodeCoordinates(elementNodes(:,2),2)-nodeCoordinates(elementNodes(:,1),2);
% % % % % % %     za=nodeCoordinates(elementNodes(:,2),3)-nodeCoordinates(elementNodes(:,1),3);
% % % % % % %    
% % % % % % %     Ldef=sqrt(xa.^2+ya.^2); %%This is the deformed length of the beam for angle computation
% % % % % % %     C=xa./Ldef; %%Global Cosine of element
% % % % % % %     S=ya./Ldef;  %%Global Sine of element
% % % % % % %     K=PlaneFrameElStiffVector(E,A,I,L,C,S,G,EI,element,elementNodes,GDof,beamType);
    LL=(xa.*xa+ ya.*ya + za.*za).^0.5;
    l=xa./LL;
    m=ya./LL;
    n=za./LL; 

    D=sqrt(l.*l+m.*m);
    Lvx=-m./D;
    Lvy=l./D;
    Lvz=0./D;
    %Lvector2=[Lvx Lvy Lvx];
    
    Lv3x=-(l.*n./D);
    Lv3y=-(m.*n./D);
    
CXx = (xa)./LL;
CYx = (ya)./LL;
CZx = (za)./LL;
D = sqrt(CXx.*CXx + CYx.*CYx);
CXy = -CYx./D;
CYy = CXx./D;
CZy = 0;
CXz = -CXx.*CZx./D;
CYz = -CYx.*CZx./D;
CZz = D;
buildK=0
K= SpaceBeamElemStiffness_mod2(E,A,GDof,LL,Ic,G,elementNodes,element,CXx,CYx,CZx,CXy,CYy,CZy,CXz,CYz,CZz);
buildK=1
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%/ ATTRACIVE FORCES BY VAN DER WAALS FORCE %%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         sizeClose=size(closeNodes,1);
         xa=zeros(sizeClose,1); ya=zeros(sizeClose,1); LL=zeros(sizeClose,1);
         CC=zeros(sizeClose,1); SS=zeros(sizeClose,1);
         xa=nodeCoordinates(closeNodes(:,2),1)-nodeCoordinates(closeNodes(:,1),1);
         ya=nodeCoordinates(closeNodes(:,2),2)-nodeCoordinates(closeNodes(:,1),2);
         za=nodeCoordinates(closeNodes(:,2),3)-nodeCoordinates(closeNodes(:,1),3);

            ll=(xb.*xb + yb.*yb + zb.*zb).^0.5;
            CXx=xb./ll;
            CYx=yb./ll;
            CZx=zb./ll;
                
            [Ac,vdwk]=ConnectionStiffness_3D(rout,j,closeNodes,GDof,beamType,CXx,CYx,CZx);
             K=K+Ac;

    pdof=1:numberBeams;
    xdof=GDof-6*(pdof-1)-5; %% horizontal displacement DOF for each Beam
    ydof=GDof-6*(pdof-1)-4; %% vertical displacement DOF for each beam
    zdof=GDof-6*(pdof-1)-3;
    theta1dof=GDof-6*(pdof-1)-2;
    theta2dof=GDof-6*(pdof-1)-1;
    theta3dof=GDof-6*(pdof-1)-0;
    alldof=[xdof ydof zdof theta1dof theta2dof theta3dof];
    prescribedDof=[alldof]';
    entersolution=0
    displacements=solution(GDof,prescribedDof,K,force,reactionForce);
    leavesolution=1
    U=displacements+UU;
% % % % %              
% % % % %     pdof=1:numberBeams-1;
% % % % %     mdof=GDof-3*(pdof-1); %% moment DOF for each Beam
% % % % %     hdof=GDof-2-(3*(pdof-1)); %% horizontal displacement DOF for each beam
% % % % %     vdof=GDof-1-(3*(pdof-1));
% % % % %     
% % % % %     hhdof=3*(compressedNodes-1)+1; %% No horizontal motion of nodes touching compression platen
% % % % %     vvdof=3*(compressedNodes-1)+2; %% Vertical displacement due to compression
% % % % %     
% % % % %     alldof=[mdof hdof vdof hhdof' vvdof'];
% % % % %     prescribedDof=alldof';
% % % % %            
% % % % %     displacements=solutionComp(GDof,prescribedDof,K,forceCompress,UU); %%Now using the nodal forces to compute global displacements

    reactionForce=K*U;
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     
    
    for j=1:size(compressedNodes,1);
             activeForce(j)=K(3*(compressedNodes(j)-1)+2,:)*U;
    end
 
    compressiveLoad(t)=sum(activeForce);
    
    ii=1:nodeCount; %%Relcates Nodes
        nodeCoordinates(ii,1)=nodeCoordinates(ii,1)+U(3*(ii-1)+1);
        nodeCoordinates(ii,2)=nodeCoordinates(ii,2)+U(3*(ii-1)+2);
        nodeCoordinates(ii,3)=nodeCoordinates(ii,3)+U(3*(ii-1)+3);

%ContinuousPlot=1;
         
        if ContinuousPlot==1;
           % CNTPlot(fname,avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad)
        CNTPlot_Temp(fname, span, numberBeams,nodeCount,nodeCoordinates,t,title, T, elementNodes,TempMax)
        end   
    
 end %% End Compression Loops

t=t+1; 
%CNTPlot(avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad)
CNTPlot_Temp(fname, span, numberBeams,nodeCount,nodeCoordinates,t,title, T, elementNodes,TempMax)

 
 

[Fx,Fy,Fx2,Fy2,F_axial, F_axial2,F_transverse,F_transverse2] = ForceComp(K,U,element,numberBeams, nodeCoordinates,closeNodes,elementNodes);

%%AxialStress=F_axial./A;

f=1:numberBeams:element-numberBeams;
% 
% subplot(211)
%     plot(Fy(f))
% subplot(212)
%     plot(Fx(f))



%if ContinuousPlot==0;
%CNTPlot(avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title)%%Plots final CNT Morphology
%end

toc()


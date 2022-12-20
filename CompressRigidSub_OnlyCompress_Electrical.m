% E=1E12; rout=12.5E-9; ri=7.5E-9; A=pi*(rout^2-ri^2); I=pi/2*(rout^4-ri^4); EA=E*A; EI=E*I;
% G=5e11;

load('\\tsclient\MattMaschmann\Documents\MATLAB\ImageFolder\TwoSide\Fast400a_.mat')
%load('C:\users\MattMaschmann\Documents\MATLAB\ImageFolder\TwoSide\400b_.mat')
steps=t; %OVER-RIDING
title=input('Name of Run ','s')
fname = '\\tsclient\MattMaschmann\Documents\MATLAB\ImageFolder\TwoSide\';
%fname = 'C:\users\MattMaschmann\Documents\MATLAB\ImageFolder\TwoSide\';

tgrowth=t;
totalCompress=250;

Vleft=0; Vright=10; %Voltages
conductivity=3.3e8; % Electrical conductivity of 1 CNT ohm^-1 / meter
electrodeWidth=3e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              COMPRESSION OF CNT FOREST                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GDof=3*(size(nodeCoordinates,1));%% Adding new node at center of fiber
UU=zeros(GDof,1); forceCompress=zeros(GDof,1); U=zeros(GDof,1); compressiveLoad=zeros(t,1);internalForce=0;
%nodeCount=nodeCount+1; %% This will serve as the center fiber node
V=zeros(GDof/3,1); i=zeros(GDof/3,1); Resistance=0;

for e=element-numberBeams+1:element;
        L(e)=(rate(e-element+numberBeams));   A(e)=(A(e-element+numberBeams));I(e)=I(e-element+numberBeams)';
end

topNode=max(nodeCoordinates(:,2));
bottomNode=min(nodeCoordinates(:,2));
closeCounter=0; clear Resistance; clear totalcurrent;
for compress=1:totalCompress;
    U=zeros(GDof,1); forceCompress=zeros(GDof,1); UU=sparse(GDof,1,0,GDof,1); UUU=sparse(GDof,1,0,GDof,1); activeForce=zeros(t,1); totalU=0;
    K=sparse(5,5,0,GDof,GDof);  AA=sparse(5,5,0,GDof,GDof); KK=sparse(5,5,0,GDof,GDof); JJ=sparse(5,5,0,GDof,GDof); internalForce=sparse(1,1,0,1,GDof); forcetemp=zeros(GDof,1);
    current=zeros(nodeCount,1); compforce=zeros(GDof,1); 
    
    compress
    t=compress+steps;
    %[crossConnect, JJ]=fiberStiffnessSpan2(nodeCount,numberBeams,GDof,nodeCoordinates)
    [crossConnect, JJ,nodeCoordinates]=fiberStiffnessSpan3(nodeCount,numberBeams,GDof,nodeCoordinates); %% JJ is fiber stiffness
    [UU,compressedNodes,compressDisp] = nodeCompression(topNode,nodeCoordinates,compress,GDof); %% Finds nodes to be compressed
    %[closeNodes]=FindCloseNodes_Range(nodeCoordinates,nodeCount,numberBeams)  ;
    %FindCloseNodesSparse6(nodeCoordinates,nodeCount,numberBeams,growthtime,rmax,span) ;  
    %[closeNodes] =FindCloseNodesSparse4(nodeCoordinates,nodeCount,numberBeams,t,rmax) ;
    %[closeNodes] =FindCloseNodesSparse5(nodeCoordinates,nodeCount,numberBeams,t,rmax,span) ;
    [closeNodes]=FindCloseNodes_Range_twoSide(nodeCoordinates,nodeCount,numberBeams)  ;
 

    xa=nodeCoordinates(elementNodes(:,2),1)-nodeCoordinates(elementNodes(:,1),1);
    ya=nodeCoordinates(elementNodes(:,2),2)-nodeCoordinates(elementNodes(:,1),2);
    
    Ldef=sqrt(xa.^2+ya.^2); %%This is the deformed length of the beam for angle computation
    C=xa./Ldef; %%Global Cosine of element
    S=ya./Ldef;  %%Global Sine of element
    K=PlaneFrameElStiffVector(E,A,I,L,C,S,G,EI,element,elementNodes,GDof,beamType);
    K_Elect=PlaneFrameElec(A,L,elementNodes,GDof,conductivity);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%/ ATTRACIVE FORCES BY VAN DER WAALS FORCE %%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         sizeClose=size(closeNodes,1);
         xa=zeros(sizeClose,1); ya=zeros(sizeClose,1); LL=zeros(sizeClose,1);
         CC=zeros(sizeClose,1); SS=zeros(sizeClose,1);
         xa=nodeCoordinates(closeNodes(:,2),1)-nodeCoordinates(closeNodes(:,1),1);
         ya=nodeCoordinates(closeNodes(:,2),2)-nodeCoordinates(closeNodes(:,1),2);
         
         LL=sqrt((xa.*xa+ya.*ya));
         CC=xa./LL; %Cosine
         SS=ya./LL;  %Sine
        
         roo=12.5E-9; rii=7.5E-9; Ac=pi*(roo^2-rii^2); Ic=pi/2*(roo^4-rii^4);
         
         [AA,vdwk]=ConnectionStiffness(E,Ac,Ic,LL,CC,SS,G,EI,sizeClose,closeNodes,GDof,beamType);
         [Aelect]=ConnectionElect(sizeClose,closeNodes,GDof); %Contact Resistance

         K=K+AA+JJ;%+KK
         K_Elect=K_Elect+Aelect;
      
        [UUU,compressedNodesBottom] = nodeCompressionRadialBottom(bottomNode,nodeCoordinates,compress,GDof,numberBeams,steps); 
   %Note: These nodes might contribute to surface resistance
        
    %%%%% SOLUTION %%%%%%%%    
     pdof=1:numberBeams;
     mdof=GDof-3-3*(pdof-1); %% moment (angular) DOF for each Beam
     hdof=GDof-2-(3*(pdof-1)); %% horizontal displacement DOF for each beam
%     vdof=GDof-1-(3*(pdof-1));
     leftBot=numberBeams/2+1; %node number of left-most bottom CNT
     rightBot=numberBeams; %Node number of right-most bottom CNT

    BotHorDof=3*(compressedNodesBottom-1)+1; %NODES AT BOTTOM SURFACE
    BotVertDof=3*(compressedNodesBottom-1)+2;
    
    hhdof=3*(compressedNodes-1)+1; %% No horizontal motion of nodes touching compression platen
    vvdof=3*(compressedNodes-1)+2; %% Vertical displacement due to compression
    
    totaldofs=[BotHorDof' BotVertDof' hhdof' mdof hdof vvdof'];
    displacedDofs=[BotHorDof' BotVertDof' hhdof' vvdof'];
  
    prescribedDofs=totaldofs'; 
    displacedNodes=[compressedNodes; compressedNodesBottom];
  
    activeDof=setdiff([1:GDof]',prescribedDofs);
    displacements=full(UU)+full(UUU);
  
     %U(activeDof)=K(activeDof,activeDof)\-forcetemp(activeDof);
     compforce=-(K(activeDof,displacedDofs)*displacements(displacedDofs));
     ComputeU=0
     tic()
     U(activeDof)=K(activeDof,activeDof)\compforce;
     toc()
     ComputeU=1
    %%U(freenodes)=K(freenodes,freenodes)\-(K(freenodes,13:14)*U(13:14)');
    
    disp=[U]+[displacements];
    force(displacedDofs)=K(displacedDofs,:)*disp;
    
    compressForce=K*UU;
    
      
 %%%%%%%%%%%%%%%%%%%%%%% Electrical SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %i=zeros(GDof/3,1);
     %m=min(nodeCoordinates(:,2)); %Nodes contacting bottom Surf
     i_touch=find(nodeCoordinates(:,2) < (min(nodeCoordinates(:,2))+2*10^-9))
     contacts=size(i_touch,1);
     nodecountleft=0; nodecountright=0; electContact=0;% leftContacts(1)=numberBeams/2+1; rightContacts(1)=numberBeams-5;
     for i=1:contacts
         if nodeCoordinates(i_touch(i),1)<electrodeWidth
             nodecountleft=nodecountleft+1;
             leftContacts(nodecountleft)=i_touch(i);
         end
         if nodeCoordinates(i_touch(i),1)> (span-electrodeWidth)
             nodecountright=nodecountright+1;
             rightContacts(nodecountright)=i_touch(i);
         end
     end
         
     if nodecountleft>0 && nodecountright>0
     electContact=1;
     V=zeros(nodeCount,1);   
     V(leftContacts)=Vleft; %% Bottom nodes at 0 
     V(rightContacts)=Vright; %% Right Contacts
  
     prescribedDof_Elect=[leftContacts rightContacts]';
     [V]=solutionElect(GDof,prescribedDof_Elect,K_Elect,V,current,nodeCount,numberBeams);
     %V=V+Vactive;
     
     currents=K_Elect*V;
     totalcurrent(compress)=sum(currents(rightContacts));
     Resistance(compress)=(abs(Vright-Vleft))/totalcurrent(compress)
     %currentmax=max(current)
     end
    
    for j=1:size(compressedNodes,1);
             activeForce(j)=K(3*(compressedNodes(j)-1)+2,:)*U;
    end
 
%     compressiveLoad(t)=sum(activeForce);
%     xx=1:t;
%     plot(xx,compressiveLoad,'o')             
%     caxis([0, 10]);
%              %axis([0 axesBound*1e6 0 axesBound*1e6]);
%              figureHandle = gcf;
%              set(gca,'FontSize',22)
%              set(findall(gcf,'type','text'),'FontSize',22)
%              xlabel('time (\mum)');
%              ylabel('Force (\mum)');
%              plotname=strcat('Force',num2str(t));
%              saveas(gcf,fullfile(fname,plotname),'png');
    
    for ii=1:nodeCount; %%Relocates Nodes
        nodeCoordinates(ii,1)=nodeCoordinates(ii,1)+disp(3*(ii-1)+1);
        nodeCoordinates(ii,2)=nodeCoordinates(ii,2)+disp(3*(ii-1)+2);
    end
VMax=max(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% RESISTANCE PLOT  %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if electContact>0 && rem(compress,5)==0;
%        plot(compressDisp*(1:compress)*10^6,abs(Resistance),'o', 'Color','red')
%        axis([0 10 1 1000]);
%        xlabel('Compression (\mum)');
%        ylabel('Resistance (Ohms)');
%        name=char('Resistance');
%        plotname=strcat(title,name,num2str(t));
%        hold on
%       % saveas(gcf,fullfile(fname,plotname),'png');
%         %close()
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        if electContact>0 && rem(compress,50)==0;
           %CNTPlotFastCompress(fname,avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,tgrowth)
           CNTPlot_Elect(fname, span, numberBeams,nodeCount,nodeCoordinates,t,title, V, elementNodes,VMax,compressDisp,compress,Resistance) 
           %CNTPlotCompress(fname,avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,bottomNode)
        end   
        
        if rem(compress,20)==0
            save([fullfile(fname,titlename),'.mat']); % Saves workspace for further manipulation
        end
    
 end %% End Compression Loops


 



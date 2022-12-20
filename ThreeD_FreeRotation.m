clear
E=1e12; rout=5e-9; ri=3.5e-9; A=pi*(rout^2-ri^2); Ic=pi/2*(rout^4-ri^4); EA=E*A; EI=E*Ic; G=5e11; 
title=input('Name of Run ','s');
fname = 'C:\Users\MattMaschmann\Documents\MATLAB\ImageFolder\Junk';
%fname ='/home/mmaschma/data/3Dcode/Output/test10';
steps=1;% Number of time steps
h_span=2;% Span of the substrate lengthwise(meters)
v_span=2;% Span of the substrate breadthwise(meters)
numberBeamsx=10;% Number of CNTs to be modeled lengthwise
numberBeamsy=10;% Number of CNTs to be modeled breadthwise

truncatedTime=50;

ang_stdev=10;% Choose 0,5,10, percent for parametric study
rate_stdev=10;% Choose 0, 10, 20 percent for parametric study
avgRate=60e-9;% Average growth per time step (meters)
numberBeams=numberBeamsx*numberBeamsy;% Total number of Beams on the Substrate
ContinuousPlot=1; plotint=1;% 0=plotting off.  1=plotting on
PeriodicBoundary=0;% 0=off, 1=on
beamType=0;% 0=Euler Beam ; 1=Timoshenko Beam
totalCompress=0;
compressiveLoad=0;
element=numberBeams;
nodeCount=2*numberBeams;
contactingNodes=zeros(steps,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%/ nucleating CNTs with distributed properties according to inputs%%%%/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,ro,rinn,phi]=nucleate3D_MMm(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rout,ri,rate_stdev,ang_stdev);

A=pi*(ro.^2-ri.^2)';
Ic=pi/2*(ro.^4-ri.^4)'; 
EA=E.*A; EI=E.*Ic;
L=rate'; plotcount=0;
reactionForce=zeros(6*(numberBeams*2),1);% Element force vector
for t=1:steps %% Growing CNT forest for a quantity of "steps" time steps
    if t>1
        t
        for e=element-numberBeams+1:element;
        e;   
        L(e)=(rate(e-element+numberBeams));
        A(e)=(A(e-element+numberBeams));
        Ic(e)=Ic(e-element+numberBeams)';
        end
    end

    GDof=6*nodeCount;  
    %%closeNodes=FindCNTMeshGrid(nodeCoordinates,nodeCount,numberBeams,t,GDof) 
     rmax=max(ro);
     entercloseNodes=0
     [gap, sep,closeNodes] =FindCloseNodesSparse2(nodeCoordinates,nodeCount,numberBeamsx,numberBeamsy,t,rmax)  ; 
    %[gap, sep,closeNodes] =FindCloseNodes3DTruncated(nodeCoordinates,nodeCount,numberBeamsx,numberBeamsy,t,steps,oldCloseNodes,truncatedTime); 
     %[gap, sep,closeNodes] =FindCloseNodesSparse(nodeCoordinates,nodeCount,numberBeams,t,rmax)   ; 
     size(closeNodes)
     leaveCloseNodes=1
    %%closeNodes=FindCloseNodesEachCNT(nodeCoordinates,nodeCount,numberBeams,t,GDof)  
      
    U=zeros(GDof,1) ;
    force=zeros(GDof,1); %%initializing variables
%     K=zeros(GDof,GDof);
%     Ab=zeros(GDof,GDof);
    K=sparse(5,5,0,GDof,GDof);
    Ab=sparse(5,5,0,GDof,GDof);
%    Ac=sparse(5,5,0,GDof,GDof);
    AA=sparse(1,1,0,6,6);
    i=0 ;
    ii=zeros(144,1);
    iii=zeros(144,1);
    j=0 ;
    jj=zeros(144,1);
    jjj=zeros(144,1);
    xa=zeros(element,1);
    ya=zeros(element,1);
    za=zeros(element,1);
    C=zeros(element,1);
    s=zeros(element,1);

        
    xa=nodeCoordinates(elementNodes(:,2),1)-nodeCoordinates(elementNodes(:,1),1); % x component of the length of the element where elementNodes is the connectivity matrix   
    ya=nodeCoordinates(elementNodes(:,2),2)-nodeCoordinates(elementNodes(:,1),2); % y component of the length of the element 
    za=nodeCoordinates(elementNodes(:,2),3)-nodeCoordinates(elementNodes(:,1),3); % z componet of the length of the element
    
    if PeriodicBoundary>0;
        %for h=1:element;
           if abs(xa)>h_span/2;
                if nodeCoordinates(elementNodes(h,1),1)>nodeCoordinates(elementNodes(h,2),1) ;
                    xa=h_span+nodeCoordinates(elementNodes(h,2),1)-nodeCoordinates(elementNodes(h,1),1);
                    else
                        xa=h_span-nodeCoordinates(elementNodes(h,2),1)+nodeCoordinates(elementNodes(h,1),1);
                    end
            end
        %end
    end
      
      NDF=6; % number of DOFs per node
      NTT=nodeCount*NDF;% total number of DOFs
      %K=zeros(NTT,NTT);
      tt=zeros(3,3);
      localx=zeros(1,3);
      localy=zeros(1,3);
      localz=zeros(1,3);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%/ ATTRACIVE FORCES BY VAN DER WAALS ATTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%/    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
    

%if size(closeNodes,1)==0 
%   closeNodes=0;%%[nodeCount,0] closeNodes==[nodeCount,0]
%else
  
if closeNodes ~= 0  
    
    NumberBars=size(closeNodes,1);
    xb=nodeCoordinates(closeNodes(:,2),1)-nodeCoordinates(closeNodes(:,1),1);
    yb=nodeCoordinates(closeNodes(:,2),2)-nodeCoordinates(closeNodes(:,1),2);
    zb=nodeCoordinates(closeNodes(:,2),3)-nodeCoordinates(closeNodes(:,1),3);
    ndf=3;%number of DOFs per node in the bar element
    ntt=2*NumberBars*ndf;%total number of DOFs

    if PeriodicBoundary>0;
                %for hh=1:NumberBars;
                   if abs(xb)>h_span/2;
                        if nodeCoordinates(closeNodes(hh,1),1)>nodeCoordinates(closeNodes(hh,2),1); 
                            xb=h_span+nodeCoordinates(closeNodes(hh,2),1)-nodeCoordinates(closeNodes(h,1),1);
                            else
                                xb=h_span-nodeCoordinates(closeNodes(hh,2),1)+nodeCoordinates(closeNodes(hh,1),1);
                            end
                    end
                %end
            end
    
         Ab=sparse(1,1,0,GDof,GDof);
            ll=(xb.*xb + yb.*yb + zb.*zb).^0.5;
            CXx=xb./ll;
            CYx=yb./ll;
            CZx=zb./ll;
                
[Ac,vdwk]=ConnectionStiffness_3D(rout,j,closeNodes,GDof,beamType,CXx,CYx,CZx);
    K=K+Ac;        
end
    contactingNodes(t)=size(closeNodes,1);
    
% % bottomNodes=nodeCount-numberBeams+1:nodeCount;
% % rotation=RotationStiffness(element,bottomNodes,GDof);
% %  
% % K=K+rotation;


for ff=1:numberBeams
    ff;
    %forces(GDof-12*numberBeams+6*(ff-1)+1)=E*A(ff)*sin(ang(ff))*cos(phi(ff));
    ii(ff)=GDof-12*numberBeams+6*(ff-1)+1;
    jj(ff)=1;
    kkforce(ff)=E*A(ff)*sin(ang(ff))*cos(phi(ff));
    
    %forces(GDof-12*numberBeams+6*(ff-1)+2)=E*A(ff)*sin(ang(ff))*sin(phi(ff));
    ii(numberBeams+ff)=GDof-12*numberBeams+6*(ff-1)+2;
    jj(numberBeams+ff)=1;
    kkforce(numberBeams+ff)=E*A(ff)*sin(ang(ff))*sin(phi(ff));
    
    %forces(GDof-12*numberBeams+6*(ff-1)+3)=E*A(ff)*cos(ang(ff));
    ii(2*numberBeams+ff)=GDof-12*numberBeams+6*(ff-1)+3;
    jj(2*numberBeams+ff)=1;
    kkforce(2*numberBeams+ff)=E*A(ff)*cos(ang(ff));
    
end
    force=sparse(ii',jj',kkforce',GDof,1);
    %force=sparse(forces);

    pdof=1:numberBeams;
    xdof=GDof-6*(pdof-1)-5; %% horizontal displacement DOF for each Beam
    ydof=GDof-6*(pdof-1)-4; %% vertical displacement DOF for each beam
    zdof=GDof-6*(pdof-1)-3;
    theta1dof=GDof-6*(pdof-1)-2;
    theta2dof=GDof-6*(pdof-1)-1;
    theta3dof=GDof-6*(pdof-1)-0;
    alldof=[xdof ydof zdof];% theta1dof theta2dof theta3dof];
    prescribedDof=[alldof]';
    entersolution=0
    U=solution(GDof,prescribedDof,K,force,reactionForce);
    leavesolution=1
    %U=displacements;
    deformed_coord=zeros(nodeCount,2);
    
    %reactionForce=zeros(GDof+numberBeams*6, 1);
  %      activeDof=setdiff([1:GDof]',[prescribedDof]);
    activeDof=setdiff([1:GDof]',[prescribedDof]);
    reactionForce=K(activeDof,activeDof)*U(activeDof);
    
    for i=size(reactionForce,1)+1:size(reactionForce,1)+6*numberBeams
        reactionForce(i)=reactionForce(i-numberBeams);
    end
    
    %%%% REPOSITIONING NODES %%%%/
    
    %%for ii=1:nodeCount-numberBeams %%Shifts Upper Nodes Up
        ii=1:nodeCount-numberBeams; %%Shifts Upper Nodes Up
        nodeCoordinates(ii,1)=nodeCoordinates(ii,1)+U(6*(ii-1)+1);
        nodeCoordinates(ii,2)=nodeCoordinates(ii,2)+U(6*(ii-1)+2);
        nodeCoordinates(ii,3)=nodeCoordinates(ii,3)+U(6*(ii-1)+3);
       
        if PeriodicBoundary>0;
            for ii=1:nodeCount-numberBeams;       
                if nodeCoordinates(ii,1)<0;
                    nodeCoordinates(ii,1)=nodeCoordinates(ii,1)+h_span;
                end
                 if nodeCoordinates(ii,1)>h_span;
                    nodeCoordinates(ii,1)=nodeCoordinates(ii,1)-h_span;
                 end
            end
        end
        
  
        jj=nodeCount-numberBeams+1:nodeCount;%% shifts surface nodes up
        nodeCoordinates(jj,1)=growthNode(jj-nodeCount+numberBeams,1);
        nodeCoordinates(jj,2)=growthNode(jj-nodeCount+numberBeams,2);
        nodeCoordinates(jj,3)=growthNode(jj-nodeCount+numberBeams,3);
        elementNodes(element,1)=element;
        elementNodes(element,2)=nodeCount;
    %%end
    
    currentCount=nodeCount;
    
     for kk=currentCount+1:currentCount+numberBeams
%%        if t<steps ;
             nodeCount=nodeCount+1;
             element=element+1 ;
 %%        end
         %%kk=currentCount+1:currentCount+numberBeams
         nodeCoordinates(kk,1)=nucleationSite(kk-currentCount,1);
         nodeCoordinates(kk,2)=nucleationSite(kk-currentCount,2);
         nodeCoordinates(kk,3)=nucleationSite(kk-currentCount,3);
         elementNodes(element,1)=element;
         elementNodes(element,2)=nodeCount;
    end
    
%reactionForce = ForceComp(K,U,element,numberBeams, nodeCoordinates,closeNodes,elementNodes); %To be used for internal residual stress
%reactionForce=sparse(reactionForce);
if ContinuousPlot==1
   plotcount=plotcount+1;
   if plotcount==plotint;
   CNTPlotBW(fname,avgRate, steps, h_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,numberBeamsx,numberBeamsy);
   plotcount=0;
   end
end
    
%save([fullfile(fname,title),'.mat']) % Saves workspace for further manipulation

end
f=1:numberBeams:element-numberBeams;
save([fullfile(fname,title),'.mat']) % Saves workspace for further manipulation
if ContinuousPlot==0;
CNTPlot(fname, avgRate, steps, h_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad);%%Plots final CNT Morphology
end




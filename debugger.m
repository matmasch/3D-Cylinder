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
ang_stdev=5;% Choose 0,5,10, percent for parametric study
rate_stdev=7;% Choose 0, 10, 20 percent for parametric study
avgRate=60e-9;% Average growth per time step (meters)
numberBeams=numberBeamsx*numberBeamsy;% Total number of Beams on the Substrate
ContinuousPlot=1;% 0=plotting off.  1=plotting on
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
[elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,ro,rinn,phi]=nucleate3D_MMm(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rout,ri,rate_stdev,ang_stdev)

A=pi*(ro.^2-ri.^2)';
Ic=pi/2*(ro.^4-ri.^4)'; 
EA=E.*A; EI=E.*Ic;
L=rate';
reactionForce=zeros(6*numberBeams,1);% Element force vector
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
     
     [gap, sep,closeNodes]=FindCloseNodesSparse(nodeCoordinates,nodeCount,numberBeams,t,rmax); 
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
   
%for il=1:element;
  
   
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
    
    %Lvector3=[Lv3x Lv3y D];
    %localx=[l m n];
    
    %localy=Lvector2;
    
    %localz=Lvector3;
    
    
    k1 = E.*A./LL;
    k2 = 12.*E.*Ic./(LL.*LL.*LL);
    k3 = 6*E.*Ic./(LL.*LL);
    k4 = 4*E.*Ic./LL;
    k5 = 2*E.*Ic./LL;
    k6 = 12*E.*Ic./(LL.*LL.*LL);
    k7 = 6*E.*Ic./(LL.*LL);
    k8 = 4*E.*Ic./LL;
    k9 = 2*E.*Ic./LL;
    k10 = G.*Ic./LL;
    
    kk=sparse( 1, 1, 0, GDof, GDof)
 for il=1:element;  
    ii=zeros(144,1);
    jj=zeros(144,1);
    [Lambda]=[l(il) m(il) n(il); Lvx(il) Lvy(il) Lvz(il); Lv3x(il) Lv3y(il) D(il)];
    T=zeros(12,12); % Transformation matrix
    R = [Lambda zeros(3,9); zeros(3) Lambda zeros(3,6);zeros(3,6) Lambda zeros(3);zeros(3,9) Lambda];
     

    a=[k1(il) 0 0; 0 k2(il) 0; 0 0 k6(il)];
    b=[ 0 0 0;0 0 k3(il); 0 -k7(il) 0];
    c=[k10(il) 0 0;0 k8(il) 0; 0 0 k4(il)];
    d=[-k10(il) 0 0;0 k9(il) 0;0 0 k5(il)];
    k = [a b -a b;b' c b d; (-a)' b' a -b;b' d' (-b)' c];

    kj=R'*k*R; %Transformed stiffness matrix
    for ppp=1:144;
        ksp(ppp)=kj(ppp);
    end
    ii=[(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF], (elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF],(elementNodes(il,1)-1)*NDF+[1:NDF] , (elementNodes(il,2)-1)*NDF+[1:NDF]];
    oo=(elementNodes(il,1)-1)*NDF; ww=(elementNodes(il,2)-1)*NDF;
    
    jj=[oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+4 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+5 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 oo+6 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+4 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+5 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6 ww+6];
    kk=sparse( jj, ii, ksp, GDof, GDof);
    
    K=K+kk; % assemble the global stiffness matrix
end


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
Kez = SpaceBeamElemStiffness(E,A,GDof,LL,Ic,G,elementNodes,element,CXx,CYx,CZx,CXy,CYy,CZy,CXz,CYz,CZz);
NumberBars=size(closeNodes,1);

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
                        if nodeCoordinates(closeNodes(hh,1),1)>nodeCoordinates(closeNodes(hh,2),1) 
                            xb=h_span+nodeCoordinates(closeNodes(hh,2),1)-nodeCoordinates(closeNodes(h,1),1);
                            else
                                xb=h_span-nodeCoordinates(closeNodes(hh,2),1)+nodeCoordinates(closeNodes(hh,1),1);
                            end
                    end
                %end
            end
    
    closenodesflag=1
        size(closeNodes)
        Ab=sparse(1,1,0,GDof,GDof);
        Abadd=sparse( jj, ii, ksp, GDof, GDof);
         
            %lvector1=[xb yb zb(];
            ll=(xb.*xb + yb.*yb + zb.*zb).^0.5;
            CXx=xb./ll;
            CYx=yb./ll;
            CZx=zb./ll;
            
            for j=1:NumberBars;
                oo=0; ww=0;
            T = [CXx(j)*CXx(j) CXx(j)*CYx(j) CXx(j)*CZx(j) ; CYx(j)*CXx(j) CYx(j)*CYx(j) CYx(j)*CZx(j) ; CZx(j)*CXx(j) CZx(j)*CYx(j) CZx(j)*CZx(j)];

            BB=10000*[T -T ; -T T];
            %BB=r'*AA*r;
            
            iii=[(closeNodes(j,1)-1)*6+[1:ndf],(closeNodes(j,2)-1)*6+[1:ndf],(closeNodes(j,1)-1)*6+[1:ndf],(closeNodes(j,2)-1)*6+[1:ndf],(closeNodes(j,1)-1)*6+[1:ndf],(closeNodes(j,2)-1)*6+[1:ndf],(closeNodes(j,1)-1)*6+[1:ndf],(closeNodes(j,2)-1)*6+[1:ndf],(closeNodes(j,1)-1)*6+[1:ndf],(closeNodes(j,2)-1)*6+[1:ndf],(closeNodes(j,1)-1)*6+[1:ndf],(closeNodes(j,2)-1)*6+[1:ndf]];
            oo=(closeNodes(j,1)-1)*6; ww=(closeNodes(j,2)-1)*6;
            jjj=[oo+1 oo+1 oo+1 oo+1 oo+1 oo+1 oo+2 oo+2 oo+2 oo+2 oo+2 oo+2 oo+3 oo+3 oo+3 oo+3 oo+3 oo+3 ww+1 ww+1 ww+1 ww+1 ww+1 ww+1 ww+2 ww+2 ww+2 ww+2 ww+2 ww+2 ww+3 ww+3 ww+3 ww+3 ww+3 ww+3];
            
                for ind=1:36;
                   BBB(ind)=BB(ind);
                end
            Abadd=sparse( jjj, iii, BBB, GDof, GDof);
            Ab=Ab+Abadd;
        end
end
    contactingNodes(t)=size(closeNodes,1);
    K=K+Ab;
    
    Ab=0;
    Abadd=0;
 
        
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
    alldof=[xdof ydof zdof theta1dof theta2dof theta3dof];
    prescribedDof=[alldof]';
    
    displacements=solution(GDof,prescribedDof,K,force,reactionForce);
    U=displacements;
    deformed_coord=zeros(nodeCount,2);
    
    reactionForce=zeros(GDof, GDof);
    
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
    
reactionForce = ForceComp(K,U,element,numberBeams, nodeCoordinates,closeNodes,elementNodes); %To be used for internal residual stress
reactionForce=sparse(reactionForce);
if ContinuousPlot==1
   CNTPlot(fname,avgRate, steps, h_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad);
end
    
    
 


save([fullfile(fname,title),'.mat']) % Saves workspace for further manipulation

%[Fx,Fy,Fx2,Fy2,F_axial, F_axial2,F_transverse,F_transverse2,reactionForce] = ForceComp(K,U,element,numberBeams, nodeCoordinates,closeNodes,elementNodes);

%%AxialStress=F_axial./A;
end
f=1:numberBeams:element-numberBeams;
% 
% subplot(211)
%     plot(Fy(f))
% subplot(212)
%     plot(Fx(f))



if ContinuousPlot==0;
CNTPlot(fname, avgRate, steps, h_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad);%%Plots final CNT Morphology
end




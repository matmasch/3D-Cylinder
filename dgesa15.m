% To show 2D built-up structures
% Case: a 2D 16-bar truss
clear all, close all

% (I) Pre-Processor
dx=10; dy=10; % ft
XYZ=[0 0 0;0 dy 0;0 2*dy 0;0 3*dy 0;
    dx 0 0;dx dy 0;dx 2*dy 0;dx 3*dy 0;
    1.65*dx 2.6*dy 0;2.3*dx 3*dy 0;
    2.95*dx 2.6*dy 0;3.6*dx 3*dy 0; 
    4.25*dx 2.6*dy 0;4.9*dx 3*dy 0;
    4.9*dx 2*dy 0;4.9*dx dy 0;4.9*dx 0 0;
    5.9*dx 3*dy 0;5.9*dx 2*dy 0;5.9*dx dy 0;5.9*dx 0 0] %nodal coordinates
CONNEC=[1 2;1 5;2 3;2 5;2 6;3 4;3 6;3 7;
        4 7;4 8;5 6;6 7;7 8;7 9;8 9;8 10
        9 10;9 11;10 11;10 12;11 12;11 13
        12 13;12 14;13 14;13 15;14 15;14 18;
        15 16;15 18;15 19;16 17;16 19;16 20
        17 20;17 21;18 19;19 20;20 21] % connectivity matrix
NN=size(XYZ,1) % number of nodes
NE=size(CONNEC,1)% number of elements
e=(11.6*10e6) % Young's modulus, lbf/ft^2
area=1.2/144 % cross-sectional area, ft^2
tau=0;  % pre-tension, lbf
BCC=zeros(NN,3)
BCC((1:NN),3)=1
BCC([1 5],[1 2])=1 
BCC([17 21],[1 2])=1%specify GBCs
LOAD=[11  2  -36000]  

%(II) Processor
NDF=3 % number of DOFs per node
NTT=NN*NDF % total number of DOFs
K=zeros(NTT,NTT)
F=zeros(NTT,1)
TT=zeros(3,3) % coordinate transformation
for i1=1:NE
   i1
  dxy=XYZ(CONNEC(i1,2),:)-XYZ(CONNEC(i1,1),:)  
  le=norm(dxy)
  TT(1,:)=dxy/le
  TT(2,:)=[0 0 1]
  TT(3,:)=cross(TT(1,:),TT(2,:)) 
  x1=e*area/le
  x2=tau/le
  ke=[x1 0 0 -x1 0 0;    0 x2 0 0 -x2 0;    0 0 x2 0 0 -x2;
          -x1 0 0 x1 0 0;    0 -x2 0 0 x2 0;   0 0 -x2 0 0 x2]
      
  G=zeros(6,6)
  xv=[1 2 3]
  G(xv,xv)=TT
  G(xv+3,xv+3)=TT
  ke=G'*ke*G;
  xv=[(CONNEC(i1,1)-1)*NDF+[1:NDF], (CONNEC(i1,2)-1)*NDF+[1:NDF]]
  K(xv,xv)=K(xv,xv)+ke % assemble the global stiffness matrix
end
F((LOAD(:,1)-1)*NDF+LOAD(:,2))=LOAD(:,3) %global force vector
[xi,xj]=find(BCC==1)
xv=(xi-1)*NDF+xj % location vector of fixed DOFs 
K(xv,:)=[ ]
K(:,xv)=[ ]
F(xv)=[ ] %apply GBCs

u=zeros(NTT,1)
xx=[1:NTT]' 
xx(xv)=[] 
u(xx)=K\F % nodal displacements of free DOFs 
clear xx

%(III) Post-Processor
for i1=1:NN,  XYZ1(i1,1:3)=XYZ(i1,1:3)+u((i1-1)*3+[1 2 3])';  end %compute deformed nodal coordinates

figure(1)
hold ('on'), grid('on'), view(0,90), axis([-10 70 -10 40 0 40]), xlabel('X'), ylabel('Y'), zlabel('Z')
for i1=1:NE
  xv=[CONNEC(i1,1) CONNEC(i1,2)];
  plot3(XYZ(xv,1),XYZ(xv,2),XYZ(xv,3),'.g-',XYZ1(xv,1),XYZ1(xv,2),XYZ1(xv,3),'.r-') % undeformed & deformed configurations
end
for i1=1:NN,  text(XYZ(i1,1),XYZ(i1,2),XYZ(i1,3),num2str(i1)); end %put node numbers
   
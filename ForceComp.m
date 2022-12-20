function reactionForce = ForceComp(K,U,element,numberBeams, nodeCoordinates,closeNodes,elementNodes)
%function [Fx,Fy,Fx2,Fy2,F_axial, F_axial2,F_transverse,F_transverse2,reactionForce] = ForceComp(K,U,element,numberBeams, nodeCoordinates,closeNodes,elementNodes)

reactionForce=K*U;

% xglobal=nodeCoordinates(elementNodes(:,2),1)-nodeCoordinates(elementNodes(:,1),1);
% yglobal=nodeCoordinates(elementNodes(:,2),2)-nodeCoordinates(elementNodes(:,1),2);
% 
%     Ldef=sqrt(xglobal.^2+yglobal.^2); %%This is the deformed length of the beam 
%     Cglobal=xglobal./Ldef; %%Global Cosine of element
%     Sglobal=yglobal./Ldef;  %%Global Sine of element
%     transf=[Cglobal Sglobal; -Sglobal Cglobal]; %%Transform from global to element coordinates
%     
%     
% for e=1:element-numberBeams;
%     indice=elementNodes(e,:);
%     %%[FAxial,FTransv]=[globalForce(3*(indice(1)-1)+1) globalForce(3*(indice(1)-1)+2)]
%     Fx(e)=globalForce(3*(indice(1)-1)+1); %% Global force in X-direction at node 1
%     Fx2(e)=globalForce(3*(indice(2)-1)+1); %% Global force in X-directation at node 2
%     Fy(e)=globalForce(3*(indice(1)-1)+2);
%     Fy2(e)=globalForce(3*(indice(2)-1)+2);
%     %%forceTransv_2(indice(2))=globalForce(3*(indice(2)-1)+2);
%     F_axial(e)=Fx(e)*Cglobal(e)+Fy(e)*Sglobal(e); %% Transforms global force to local force
%     F_axial2(e)=Fx2(e)*Cglobal(e)+Fy2(e)*Sglobal(e); %% Transfors global force to local force
%     
%     F_transverse(e)=-Fx(e)*Sglobal(e)+Fy(e)*Cglobal(e);
%     F_transverse2(e)=-Fx2(e)*Sglobal(e)+Fy2(e)*Cglobal(e);
% 
%     %%forcemom_1(indice(1))=globalForce(3*(indice(1)-1)+3);
%     %%forcemom_2(indice(2))=globalForce(3*(indice(2)-1)+3);
% end
%     %%F_axial=Fx_transf.*Cglobal+Fy_transf.*Sglobal;
%     
%     %%F_transverse=-Fx_transf.*Sglobal+Fy_transf.*Cglobal;







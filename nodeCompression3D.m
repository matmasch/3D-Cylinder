function [UU,compressedNodes] = nodeCompression3D(topNode,nodeCoordinates,compress,GDof)      
UU=zeros(GDof,1);
compressDisp=25e-9;

compressedNodes=find(nodeCoordinates(:,3)> topNode-compress*compressDisp);
% ii=zeros(GDof,1); jj=zeros(GDof,1); kk=zeros(GDof,1);

for i=1:size(compressedNodes,1);
    ii(i)=(6*(compressedNodes(i)-1)+3);
    %ii(i)=compressedNodes(i);
    %iilat(i)=(3*(compressedNodes(i)-1)+1);
    jj(i)=1;
    %jjlat(i)=2;
    kk(i)=topNode-compress*compressDisp-nodeCoordinates(compressedNodes(i),3);
    %kklat(i)=compressDisp;
end

UU=sparse(ii,jj,kk,GDof,1);

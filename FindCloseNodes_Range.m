function [closeNodes]=FindCloseNodes_Range(nodeCoordinates,nodeCount)  

%%% 3D find close nodes
voxel=zeros(nodeCount,2);

gap=50e-9;gapSquared=gap^2;
gridsize=1/(gap); %setting voxel size
ii=0; jj=0;
maxx=max(ceil(gridsize*nodeCoordinates(:,1))); maxy=max(ceil(gridsize*nodeCoordinates(:,2))); maxxy=maxx*maxy;

closeNodes=zeros(nodeCount,2);

counter=1;
nodenum=1:nodeCount;

[voxel]=[nodenum(:) (ceil(gridsize*nodeCoordinates(:,2))-1)*maxx  + ceil(gridsize*nodeCoordinates(:,1))  + maxxy*ceil(gridsize*nodeCoordinates(:,3))];

%[voxel]=[nodenum(:) ceil(gridsize*nodeCoordinates(:,1))+numberBeamsy*ceil(gridsize*nodeCoordinates(:,2))];

vox=voxel(:,2);
uniquevox=unique(vox);
[bincounts,vox_bin]=histc(vox, uniquevox); %vox_bin represents the bin that the voxel went into
nlist = find(bincounts > 1); % nlist is a list of nn that are repeated
occupiedvox=uniquevox(nlist);
nodelist=ismember(voxel(:,2),occupiedvox);
nodes=voxel(nodelist,:);
sortNodes=sortrows(nodes,2); %%Node Voxel
sizenodes=size(sortNodes,1);

          for s=1:sizenodes-1
              
              ss=s+1;%for ss = s+1: min(s+5,sizenodes)
                                if sortNodes(s,2)==sortNodes(ss,2)
                                 jj(counter)=sortNodes(ss,1);
                                 ii(counter)=sortNodes(s,1);
                                 dist(counter)=(   (nodeCoordinates(sortNodes(s),1)-nodeCoordinates(sortNodes(ss),1)).^2 + (nodeCoordinates(sortNodes(s),2)-nodeCoordinates(sortNodes(ss),2) ).^2 + + (nodeCoordinates(sortNodes(s),3)-nodeCoordinates(sortNodes(ss),3) ).^2);
                                    if dist(counter)<gapSquared
                                        counter=counter+1;
                                    end
                                end
                        
                 %end
          end
              
    closeNodes=[ii(:),jj(:)];
    closeNodes=unique(closeNodes,'rows');


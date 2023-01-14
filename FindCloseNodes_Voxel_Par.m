function [closeNodes]=FindCloseNodes_Voxel_Par(nodeCoordinates,nodeCount)  

%%% 3D find close nodes
voxel=zeros(nodeCount,2);

gap=30e-9;gapSquared=gap^2;
gridsize=1/(6*gap); %setting voxel size

%gridsize=1/(2*gap); %setting voxel size

ii=0; jj=0;
min_x=min(nodeCoordinates(:,1)); max_x=max(nodeCoordinates(:,1)); diffx = max_x - min_x;
min_y=min(nodeCoordinates(:,2)); max_y=max(nodeCoordinates(:,2)); diffy = max_y - min_y;
min_z=min(nodeCoordinates(:,3)); max_z=max(nodeCoordinates(:,3)); diffz = max_z - min_z;

maxx=max(ceil(gridsize*diffx));                  %%% Number of voxels in x
maxy=max(ceil(gridsize*diffy));                  %%% Number of voxels in y
maxz=max(ceil(gridsize*diffz));                  %%% Number of voxels in z

maxxy=maxx*maxy;                                %%% Total number of voxels

closeNodes=zeros(nodeCount,2);

counter=1;
nodenum=1:nodeCount;

%[voxel]=[nodenum(:) (ceil(gridsize*nodeCoordinates(:,2)-min_y))*maxx  + ceil(gridsize*(nodeCoordinates(:,1)-min_x+1e-20))  ];

[voxel]=[nodenum(:) (ceil(gridsize*(nodeCoordinates(:,3)-min_z))*maxx*maxy + ceil(gridsize*(nodeCoordinates(:,2)-min_y))*maxx  + ceil(gridsize*(nodeCoordinates(:,1)-min_x+1e-20)) )];
        %%% node    %%%Number of voxels below the current row     

%[voxel]=[nodenum(:) ceil(gridsize*nodeCoordinates(:,1))+numberBeamsy*ceil(gridsize*nodeCoordinates(:,2))];

vox=voxel(:,2);         %%% List of voxels listed per node
uniquevox=unique(vox);  %%% List of unique voxels - none repeated
[bincounts,vox_bin]=histcounts(vox, uniquevox); %vox_bin represents the bin that the voxel went into
nlist = find(bincounts > 1);    %% nlist is a list of nn that are repeated (more than one node per voxel)
occupiedvox=uniquevox(nlist);   %% Occupied voxels
nodelist=ismember(voxel(:,2),occupiedvox); %% Nodes which occupy voxels with other nodes
nodes=voxel(nodelist,:);        %% List of node number and voxel number

%% Try a parfor loop here

% % % % for s=1:size(occupiedvox,1) %% Loop through occupied voxels
% % % %     g = occupiedvox(s);     %% g is the voxel under evaluation
% % % %         inVoxel = nodes(:,2)== g;           %% inVoxel = List of ones and zeros
% % % %         nodesInVoxel = nodes(inVoxel,:);    %% nodesInVoxel = Nodes to check dist
% % % %         
% % % %         for p=1:size(nodesInVoxel,1)-1
% % % %             for pp = p+1:size(nodesInVoxel,1)
% % % %         
% % % %                 dist(counter)=(   (nodeCoordinates(nodesInVoxel(p),1)-nodeCoordinates(nodesInVoxel(pp),1)).^2 + (nodeCoordinates(nodesInVoxel(p),2)-nodeCoordinates(nodesInVoxel(pp),2) ).^2  + (nodeCoordinates(nodesInVoxel(p),3)-nodeCoordinates(nodesInVoxel(pp),3) ).^2);
% % % %                  if dist(counter)<gapSquared                           
% % % %                     jj(counter)=nodesInVoxel(pp,1);
% % % %                     ii(counter)=nodesInVoxel(p,1);
% % % %                     counter=counter+1;
% % % %                  end
% % % %             end
% % % %         end
% % % % end

closeNodes=[];

       parfor s=1 : size(occupiedvox,1)
                g = occupiedvox(s);     %% g is the voxel under evaluation
                inVoxel = nodes(:,2)== g;           %% inVoxel = List of ones and zeros
                nodesInVoxel = nodes(inVoxel,:);    %% nodesInVoxel = Nodes to check dist
                
                        for p=1:size(nodesInVoxel,1)-1
                            for pp = p+1:size(nodesInVoxel,1)
                                dist=(   (nodeCoordinates(nodesInVoxel(p),1)-nodeCoordinates(nodesInVoxel(pp),1)).^2 + (nodeCoordinates(nodesInVoxel(p),2)-nodeCoordinates(nodesInVoxel(pp),2) ).^2  + (nodeCoordinates(nodesInVoxel(p),3)-nodeCoordinates(nodesInVoxel(pp),3) ).^2);
                                if dist<gapSquared                           
                                closeNodes=[closeNodes; [nodesInVoxel(p) nodesInVoxel(pp)]];
                                end
                            end
                        end
       end
                
                
                
 


gridsize=gridsize/1.37;
%[voxel]=[nodenum(:) (ceil(gridsize*nodeCoordinates(:,2)))*maxx  + ceil(gridsize*(nodeCoordinates(:,1)-min_x+1e-20))  ];
        %%% node    %%%Number of voxels below the current row     

[voxel]=[nodenum(:) (ceil(gridsize*(nodeCoordinates(:,3)-min_z))*maxx*maxy + ceil(gridsize*(nodeCoordinates(:,2)-min_y))*maxx  + ceil(gridsize*(nodeCoordinates(:,1)-min_x+1e-20)) )];
 

vox=voxel(:,2);         %%% List of voxels listed per node
uniquevox=unique(vox);  %%% List of unique voxels - none repeated
[bincounts,vox_bin]=histcounts(vox, uniquevox); %vox_bin represents the bin that the voxel went into
nlist = find(bincounts > 1);    %% nlist is a list of nn that are repeated (more than one node per voxel)
occupiedvox=uniquevox(nlist);   %% Occupied voxels
nodelist=ismember(voxel(:,2),occupiedvox); %% Nodes which occupy voxels with other nodes
nodes=voxel(nodelist,:);        %% List of node number and voxel number
%sortNodes=sortrows(nodes,2);    %% List of occupied voxels [node voxel#]
%sizenodes=size(sortNodes,1);    %% Number of nodes that could have potential neighbors

       parfor s=1 : size(occupiedvox,1)
                g = occupiedvox(s);     %% g is the voxel under evaluation
                inVoxel = nodes(:,2)== g;           %% inVoxel = List of ones and zeros
                nodesInVoxel = nodes(inVoxel,:);    %% nodesInVoxel = Nodes to check dist
                
                        for p=1:size(nodesInVoxel,1)-1
                            for pp = p+1:size(nodesInVoxel,1)
                                dist=(   (nodeCoordinates(nodesInVoxel(p),1)-nodeCoordinates(nodesInVoxel(pp),1)).^2 + (nodeCoordinates(nodesInVoxel(p),2)-nodeCoordinates(nodesInVoxel(pp),2) ).^2  + (nodeCoordinates(nodesInVoxel(p),3)-nodeCoordinates(nodesInVoxel(pp),3) ).^2);
                                if dist<gapSquared                           
                                closeNodes=[closeNodes; [nodesInVoxel(p) nodesInVoxel(pp)]];
                                end
                            end
                        end
       end


% % % % % for s=1:size(occupiedvox,1) %% Loop through occupied voxels
% % % % %     g = occupiedvox(s);     %% g is the voxel under evaluation
% % % % %         inVoxel = nodes(:,2)== g;           %% inVoxel = List of ones and zeros
% % % % %         nodesInVoxel = nodes(inVoxel,:);    %% nodesInVoxel = Nodes to check dist
% % % % %         
% % % % %         for p=1:size(nodesInVoxel,1)-1
% % % % %             for pp = p:size(nodesInVoxel,1)    
% % % % %                 dist(counter)=(   (nodeCoordinates(nodesInVoxel(p),1)-nodeCoordinates(nodesInVoxel(pp),1)).^2 + (nodeCoordinates(nodesInVoxel(p),2)-nodeCoordinates(nodesInVoxel(pp),2) ).^2 + (nodeCoordinates(nodesInVoxel(p),3)-nodeCoordinates(nodesInVoxel(pp),3) ).^2);
% % % % %                  if dist(counter)<gapSquared                           
% % % % %                     jj(counter)=nodesInVoxel(pp,1);
% % % % %                     ii(counter)=nodesInVoxel(p,1);
% % % % %                     counter=counter+1;
% % % % %                  end
% % % % %             end
% % % % %         end
% % % % % end


              
%%%%%
    if size(closeNodes,1)==0
        closeNodes=[0,0];
    end
    closeNodes=unique(closeNodes,'rows');


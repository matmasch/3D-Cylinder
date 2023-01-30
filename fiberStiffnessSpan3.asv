function [crossConnect, JJ,nodeCoordinates]=fiberStiffnessSpan3(nodeCount,numberBeams,GDof,nodeCoordinates,length)
%erased unused quantities E, Ac, Ic, LL, G, EI
    
    alpha=1;
    epsilon=0;
    k_fiber=1e10;
    E=k_fiber;
    vdwk=k_fiber;

    ii_local=zeros(144,1); jj_local=zeros(144,1);
    
    crossConnect=zeros(numberBeams-1,2);
    
    perimeter=nodeCount-numberBeams:nodeCount-1; %Perimeter is a list of nodes residing on fiber surface
    counter=0;

%% Would it be better to create a new node? Nodecount + 1

    nodeCoordinates(nodeCount,1)=0; %generates a node at (0,0, l/2) at the last node
    nodeCoordinates(nodeCount,2)=0;
    nodeCoordinates(nodeCount,3)= length/2;   

%% CONNECTING ALL NODES WITH THE CENTER NODE OF THE FIBER

    for i=1:numberBeams;
            crossConnect(i,1)=(nodeCount-1)-numberBeams+i; %% Surface node of each CNT
            crossConnect(i,2)=nodeCount;                    % Center Node       
    end

    
    Ni=3*(crossConnect(:,1)-1);
    Nj=3*(crossConnect(:,2)-1);
        
    ii_local=[Ni+1 Ni+2 Ni+3 Nj+1 Nj+2 Nj+3 Ni+1 Ni+2 Ni+3 Nj+1 Nj+2 Nj+3 Ni+1 Ni+2 Ni+3 Nj+1 Nj+2 Nj+3 Ni+1 Ni+2 Ni+3 Nj+1 Nj+2 Nj+3 Ni+1 Ni+2 Ni+3 Nj+1 Nj+2 Nj+3 Ni+1 Ni+2 Ni+3 Nj+1 Nj+2 Nj+3]';
    jj_local=[Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3]';
    
    xa=nodeCoordinates(crossConnect(:,2),1)-nodeCoordinates(crossConnect(:,1),1);
    ya=nodeCoordinates(crossConnect(:,2),2)-nodeCoordinates(crossConnect(:,1),2);
    za=nodeCoordinates(crossConnect(:,2),3)-nodeCoordinates(crossConnect(:,1),3);
    
    L=sqrt(xa.^2+ya.^2+za.^2); %%This is the deformed length of the beam for angle computation

    %% Need to adjust to 3D spring

    %CC=xa./L; %%Global Cosine of element
    %SS=ya./L;  %%Global Sine of element

    CC=L./L; %%Global Cosine of element
    SS=L./L;  %%Global Sine of element
    
    w1=k_fiber.*CC.*CC;
    w2=k_fiber.*SS.*CC;
    w3=k_fiber.*SS.*SS;

    w1=k_fiber;     %Similarly stiff in all directions
    w2=w1;
    w3'w1;
    
   %Kg=zeros(36,numberBeams*(numberBeams-1)/2); 
   Kg=zeros(36,numberBeams); 
   
   
   Kg(1,:)= w1';
   Kg(2,:)= w2';
   Kg(3,:)= 0;
   
   Kg(4,:)= -w1';
   Kg(5,:)= -w2';
   Kg(6,:)= 0;
   
   Kg(7,:)= w2';
   Kg(8,:)= w3';
   Kg(9,:)= 0;
   
   Kg(10,:)=-w2';
   Kg(11,:)=-w3';
   Kg(12,:)=0;
   
   Kg(13,:)=0;
   Kg(14,:)=0;
   Kg(15,:)=0;
   Kg(16,:)=0;
   Kg(17,:)=0;
   Kg(18,:)=0;
   
   Kg(19,:)=-w1';
   Kg(20,:)=-w2';
   Kg(21,:)=0;   
   
   Kg(22,:)=w1';
   Kg(23,:)=w2';
   Kg(24,:)=0;     
   
   Kg(25,:)=-w2';
   Kg(26,:)=-w3';
   Kg(27,:)=0;
   
   Kg(28,:)=w2';
   Kg(29,:)=w3';
   Kg(30,:)=0;
   
   Kg(31,:)=0;
   Kg(32,:)=0;
   Kg(33,:)=0;
   Kg(34,:)=0;
   Kg(35,:)=0;
   Kg(36,:)=0;

   
   JJ=sparse(jj_local(:), ii_local(:), Kg(:), GDof, GDof);

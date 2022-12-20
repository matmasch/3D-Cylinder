function rotation=RotationStiffness(element,bottomNodes,GDof)

    ii_local=zeros(144,1); jj_local=zeros(144,1);
    Ni=6*(bottomNodes(:)-1);
    Nj=6*(bottomNodes(:)-1);
    
ii_local=[Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6]';
 
jj_local=[Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6]';

   Kg=zeros(144,size(bottomNodes,2));
   
   Kg(122,:)=1e5; % Rotational Spring Constant
   Kg(133,:)=1e5; % Rotational Spring Constant
   Kg(144,:)=1e5; % Rotational Spring Constant
   
   rotation=sparse(jj_local(:), ii_local(:), Kg(:), GDof, GDof);

       



function [Ac,vdwk]=ConnectionStiffness_3D(rout,j,closeNodes,GDof,beamType,CXx,CYx,CZx)
%erased unused quantities E, Ac, Ic, LL, G, EI
   vdwk=1000;
   
% %   alpha=1;
% %     if beamType==0; 
% %         epsilon=0;
% %     else
% %         epsilon=5000;%12.*E.*Ic./(G.*(Ac./alpha));
% %     end
% % 
% %     %vdwk=E;
% %     if rout==5e-9
% %         vdwk=273; %Spring stiffness of bar element between CNTs
% %     end
% % 
% %     if rout==12.5e-9
% %         vdwk=430;
% %     end

    ii_local=zeros(36,1); jj_local=zeros(36,1);
    
    Ni=6*(closeNodes(:,1)-1);
    Nj=6*(closeNodes(:,2)-1);
        
ii_local=[Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6]';
 
jj_local=[Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6]';

      Kg=zeros(144,size(closeNodes,1));
   
  
      Kg(1,:)=vdwk*(  CXx.^2  )';
      Kg(2,:)=vdwk*(  CXx.*CYx  )';
      Kg(3,:)=vdwk*(  CXx.*CZx  )';
      Kg(7,:)=vdwk*(  -CXx.^2  )';
      Kg(8,:)=vdwk*(  -CXx.*CYx  )';
      Kg(9,:)=vdwk*(  -CXx.*CZx  )';
      
      Kg(13,:)=vdwk*(  CXx.*CYx  )';
      Kg(14,:)=vdwk*(  CYx.^2  )';
      Kg(15,:)=vdwk*(  CYx.*CZx  )';
      Kg(19,:)=vdwk*(  -CXx.*CYx  )';
      Kg(20,:)=vdwk*(  -CYx.^2  )';
      Kg(21,:)=vdwk*(  -CYx.*CZx  )';
      
      Kg(25,:)=vdwk*(  CXx.*CZx  )';
      Kg(26,:)=vdwk*(  CYx.*CZx  )';
      Kg(27,:)=vdwk*(  CZx.^2  )';
      Kg(31,:)=vdwk*(  -CXx.*CZx  )';
      Kg(32,:)=vdwk*(  -CYx.*CZx  )';
      Kg(33,:)=vdwk*(  -CZx.^2  )';
      
      Kg(73,:)=vdwk*(  -CXx.^2  )';
      Kg(74,:)=vdwk*(  -CXx.*CYx  )';
      Kg(75,:)=vdwk*(  -CXx.*CZx  )';
      Kg(79,:)=vdwk*(  CXx.^2  )';
      Kg(80,:)=vdwk*(  CXx.*CYx  )';
      Kg(81,:)=vdwk*(  CXx.*CZx  )';
      
      Kg(85,:)=vdwk*(  -CXx.*CYx  )';
      Kg(86,:)=vdwk*(  -CYx.^2  )';
      Kg(87,:)=vdwk*(  -CYx.*CZx  )';
      Kg(91,:)=vdwk*(  CXx.*CYx  )';
      Kg(92,:)=vdwk*(  CYx.^2  )';
      Kg(93,:)=vdwk*(  CYx.*CZx  )';
      
      Kg(97,:)=vdwk*(  -CXx.*CZx  )';
      Kg(98,:)=vdwk*(  -CYx.*CZx  )';
      Kg(99,:)=vdwk*(  -CZx.^2  )';
      Kg(103,:)=vdwk*(  CXx.*CZx  )';
      Kg(104,:)=vdwk*(  CYx.*CZx  )';
      Kg(105,:)=vdwk*(  CZx.^2  )';
  
      
  
   Ac=sparse(jj_local(:), ii_local(:), Kg(:), GDof, GDof);
%%AA=sparse(ii_local(:), jj_local(:), Kg(:), GDof, GDof);

   


function [AA,vdwk]=ConnectionStiffness(E,Ac,Ic,LL,Cx,Cy,Cz,G,EI,sizeClose,closeNodes,GDof,beamType)
    
  alpha=1;
    if beamType==0; 
        epsilon=0;
    else
        epsilon=12.*E.*Ic./(G.*(Ac./alpha));
    end

    vdwk=0.05*E;

    ii_local=zeros(144,1); jj_local=zeros(144,1);
    Ni=6*(closeNodes(:,1)-1);
    Nj=6*(closeNodes(:,2)-1);
        
    ii_local=[Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6]';
    jj_local=[Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6]';
    

    w1=vdwk.*Cx.*Cx;
    w2=vdwk.*Cx.*Cy;
    w3=vdwk.*Cx.*Cz;
    w4=vdwk.*Cy.*Cy;
    w5=vdwk.*Cy.*Cz;
    w6=vdwk.*Cz.*Cz;
       
   Kg=zeros(144,sizeClose);
   
   Kg(1,:)= w1';
   Kg(2,:)= w2';
   Kg(3,:)= w3';
   
   Kg(4,:)= 0;
   Kg(5,:)= 0;
   Kg(6,:)= 0;
   
   Kg(7,:)=-w1';
   Kg(8,:)= -w2';
   Kg(9,:)= -w3';
   
   Kg(10,:)=0;
   Kg(11,:)=0;
   Kg(12,:)=0;
   
   Kg(13,:)=w2';
   Kg(14,:)=w4';
   Kg(15,:)=w5';

   Kg(16,:)=0;
   Kg(17,:)=0;
   Kg(18,:)=0;
   
   Kg(19,:)=-w2';
   Kg(20,:)=-w4';
   Kg(21,:)=-w5;   
   
   Kg(22,:)=0;
   Kg(23,:)=0;
   Kg(24,:)=0;     
   
   Kg(25,:)=w3';
   Kg(26,:)=w5';
   Kg(27,:)=w6';
   
   Kg(28,:)=0;
   Kg(29,:)=0;
   Kg(30,:)=0;
   
   Kg(31,:)=-w3';
   Kg(32,:)=-w5;
   Kg(33,:)=-w6;

   Kg(34,:)=0;
   Kg(35,:)=0;
   Kg(36,:)=0;
   
   Kg(37,:)= 0;
   Kg(38,:)= 0;
   Kg(39,:)= 0;
   
   Kg(40,:)= 0';
   Kg(41,:)= 0;
   Kg(42,:)= 0;
   
   Kg(43,:)= 0;
   Kg(44,:)= 0;
   Kg(45,:)= 0;
   
   Kg(46,:)=0;
   Kg(47,:)=0;
   Kg(48,:)=0;
   
   Kg(49,:)=0;
   Kg(50,:)=0;
   Kg(51,:)=0;
   Kg(52,:)=0;
   Kg(53,:)=0;
   Kg(54,:)=0;
   
   Kg(55,:)=0;
   Kg(56,:)=0;
   Kg(57,:)=0;   
   
   Kg(58,:)=0;
   Kg(59,:)=0;
   Kg(60,:)=0;     
   
   Kg(61,:)=0;
   Kg(62,:)=0;
   Kg(63,:)=0;
   
   Kg(64,:)=0;
   Kg(65,:)=0;
   Kg(66,:)=0;
   
   Kg(67,:)=0;
   Kg(68,:)=0;
   Kg(69,:)=0;
   Kg(70,:)=0;
   Kg(71,:)=0;
   Kg(72,:)=0;
   
   Kg(73,:)= -w1';
   Kg(74,:)= -w2';
   Kg(75,:)= -w3;
   
   Kg(76,:)= 0;
   Kg(77,:)= 0;
   Kg(78,:)= 0;
   
   Kg(79,:)= w1';
   Kg(80,:)= w2';
   Kg(81,:)= w3';
   
   Kg(82,:)=0;
   Kg(83,:)=0;
   Kg(84,:)=0;
   
   Kg(85,:)=-w2';
   Kg(86,:)=-w4';
   Kg(87,:)=-w5'

   Kg(88,:)=0;
   Kg(89,:)=0;
   Kg(90,:)=0;
   
   Kg(91,:)=w2';
   Kg(92,:)=w4';
   Kg(93,:)=w5';   
   
   Kg(94,:)=0;
   Kg(95,:)=0;
   Kg(96,:)=0;     
   
   Kg(97,:)=-w3';
   Kg(98,:)=-w5';
   Kg(99,:)=-w6;
   
   Kg(100,:)=0;
   Kg(101,:)=0;
   Kg(102,:)=0;
   
   Kg(103,:)=w3';
   Kg(104,:)=w5';
   Kg(105,:)=w6';

   Kg(106,:)=0;
   Kg(107,:)=0;
   Kg(108,:)=0;
   
   Kg(109,:)= 0;
   Kg(110,:)= 0;
   Kg(111,:)= 0;
   
   Kg(112,:)= 0;
   Kg(113,:)= 0;
   Kg(114,:)= 0;
   
   Kg(115,:)= 0;
   Kg(116,:)= 0;
   Kg(117,:)= 0;
   
   Kg(118,:)=0;
   Kg(119,:)=0;
   Kg(120,:)=0;
   
   Kg(121,:)=0;
   Kg(122,:)=0;
   Kg(123,:)=0;
   Kg(124,:)=0;
   Kg(125,:)=0;
   Kg(126,:)=0;
   
   Kg(127,:)=0;
   Kg(128,:)=0;
   Kg(129,:)=0;   
   
   Kg(130,:)=0;
   Kg(131,:)=0;
   Kg(132,:)=0;     
   
   Kg(133,:)=0;
   Kg(134,:)=0;
   Kg(135,:)=0;
   
   Kg(136,:)=0;
   Kg(137,:)=0;
   Kg(138,:)=0;
   
   Kg(139,:)=0;
   Kg(140,:)=0;
   Kg(141,:)=0;
   Kg(142,:)=0;
   Kg(143,:)=0;
   Kg(144,:)=0;

   
   %AA=sparse(jj_local(:), ii_local(:), Kg(:), GDof, GDof);
   
   AA=fsparse(jj_local(:), ii_local(:), Kg(:), [GDof,GDof]);


   


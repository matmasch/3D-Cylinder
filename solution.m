function displacements=solution(GDof,prescribedDof,K,force,reactionForce)
  
    activeDof=setdiff([1:GDof]',[prescribedDof]);

    %U=K(activeDof,activeDof)\(force(activeDof)-reactionForce(activeDof));
   
    %dk=decomposition(K(activeDof,activeDof),'chol','upper');
    %U=dk\force(activeDof);
   
   
    U = mldivide(K(activeDof,activeDof),(force(activeDof)));%-reactionForce(activeDof)));


    
    %%[U,flag,err,iter,res] = gmres(K(activeDof,activeDof), force(activeDof))
    %%U=pcg(K(activeDof,activeDof)\force(activeDof))
    %%U=umfpack(K(activeDof,activeDof),"\",force(activeDof))
    %%LUptr=umf_lufact(K(activeDof,activeDof));
    %%U=umf_lusolve(LUptr,force(activeDof));
    displacements=zeros(GDof,1);
    displacements(activeDof)=U;


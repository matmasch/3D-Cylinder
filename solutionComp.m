function displacements=solutionComp(GDof,prescribedDof,K,force,UU)
  
    activeDof=setdiff([1:GDof]',[prescribedDof]);
    
    U=K(activeDof,activeDof)\(force(activeDof)- K(activeDof,:)*UU)   ; %Should this be K(activeDof,:)?

    %%[U,flag,err,iter,res] = gmres(K(activeDof,activeDof), force(activeDof))
    %%U=pcg(K(activeDof,activeDof)\force(activeDof))
    %%U=umfpack(K(activeDof,activeDof),"\",force(activeDof))
    %%LUptr=umf_lufact(K(activeDof,activeDof));
    %%U=umf_lusolve(LUptr,force(activeDof));
    displacements=zeros(GDof,1);
    displacements(activeDof)=U;
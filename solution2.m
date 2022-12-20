function displacements=solution2(GDof,prescribedDof,K,force)
  
    activeDof=setdiff([1:GDof]',[prescribedDof])

    U=K(activeDof,activeDof)\(force(activeDof))

    displacements=zeros(GDof,1);
    displacements(activeDof)=U;


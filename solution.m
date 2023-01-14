function displacements=solution(GDof,prescribedDof,K,force,reactionForce)
  
    activeDof=setdiff([1:GDof]',[prescribedDof]);

    %U=K(activeDof,activeDof)\(force(activeDof)-reactionForce(activeDof));
   
    dk=decomposition(K(activeDof,activeDof),'chol','upper');
    U=dk\force(activeDof);

    displacements=zeros(GDof,1);
    displacements(activeDof)=U;


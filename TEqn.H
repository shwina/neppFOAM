  
	
    fvScalarMatrix TEqn
    (
      fvm::ddt(T)
    + fvm::div(phi, T)
    - fvm::laplacian(alphaEffFluid/eps, T)
    + (htc/(rhoFluid*CpF*eps))*T
    ==(htc/(rhoFluid*CpF*eps))*Tsol

    );

    TEqn.relax();
    TEqn.solve();

	

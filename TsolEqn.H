

  fvScalarMatrix TsolEqn
  (
	fvm::ddt(Tsol)
  -	fvm::laplacian(alphaEffSolid/(1-eps), Tsol)
  +	(htc/(rhoSolid*CpS*(1-eps)))*Tsol
  ==	(htc/(rhoSolid*CpS*(1-eps)))*T
  );
  
  TsolEqn.solve();



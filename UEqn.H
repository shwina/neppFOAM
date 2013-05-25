    // Construct the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
	fvm::ddt(eps, U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        sources(U)
	+ fvm::laplacian(nu, U)	
    );

    UEqn().relax();

    sources.constrain(UEqn());

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal

    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    {
        pZones.addResistance(UEqn());

        solve(UEqn() == -fvc::grad(p));

        trAU = 1.0/UEqn().A();
        trAU().rename("rAU");
    }

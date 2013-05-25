    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh

	);
    Info<< "Reading field Tsol (Solid Phase Temperature)\n" << endl;
    volScalarField Tsol
    (
        IOobject
        (
            "Tsol",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);
	
	autoPtr<incompressible::turbulenceModel> turbulence
	(
	    incompressible::turbulenceModel::New(U, phi, laminarTransport)
	);

	dimensionedScalar nu (laminarTransport.lookup("nu"));
	dimensionedScalar alphaEffFluid (laminarTransport.lookup("alphaEffFluid"));
	dimensionedScalar alphaEffSolid (laminarTransport.lookup("alphaEffSolid"));
	dimensionedScalar CpF (laminarTransport.lookup("CpF"));
	dimensionedScalar CpS (laminarTransport.lookup("CpS"));
	dimensionedScalar eps (laminarTransport.lookup("eps"));
	dimensionedScalar htc (laminarTransport.lookup("htc"));
	dimensionedScalar rhoFluid (laminarTransport.lookup("rhoFluid"));
	dimensionedScalar rhoSolid (laminarTransport.lookup("rhoSolid"));

// 
    IObasicSourceList sources(mesh);

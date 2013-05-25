if (pressureImplicitPorosity)
{
    U = trTU()&UEqn().H();
}
else
{
    U = trAU()*UEqn().H();
}
if (pimple.nCorrPISO() <= 1)
{
	UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf())
+ fvc::ddtPhiCorr(trAU, U, phi);

adjustPhi(phi, U, p);

while (pimple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn = (fvm::laplacian(trTU(), p) == fvc::div(phi));
    }
    else
    {
        tpEqn = (fvm::laplacian(trAU(), p) == fvc::div(phi));
    }

    tpEqn().setReference(pRefCell, pRefValue);

    //tpEqn().solve();
	tpEqn().solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
	if (pimple.finalNonOrthogonalIter())
    {
        phi -= tpEqn().flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();
sources.correct(U); 


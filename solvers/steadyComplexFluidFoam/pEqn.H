{
    volScalarField rAU(1.0/UEqn().A());
    //surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();
    UEqn.clear();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(HbyA) & mesh.Sf()
    );

    fvOptions.makeRelative(phiHbyA);

    adjustPhi(phiHbyA, U, p);

    // Store pressure for under-relaxation
    //p.storePrevIter();

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}

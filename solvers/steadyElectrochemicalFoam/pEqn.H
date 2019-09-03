{
    volScalarField rAU(1.0/UEqn().A());
    //surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    //if (pimple.nCorrPISO() <= 1)
    //{
        UEqn.clear();
    //}

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
     // + rAUf*fvc::ddtCorr(U, phi)
    );

    fvOptions.makeRelative(phiHbyA);

    adjustPhi(phiHbyA, U, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    K = 0.5*magSqr(U);
}

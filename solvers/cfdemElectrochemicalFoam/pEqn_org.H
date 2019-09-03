{
    //rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    //thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("rhorAUf", fvc::interpolate(rAU));
    volScalarField rAUvoidfraction("(voidfraction2|A(U))",rAU*voidfraction);
    //surfaceScalarField voidfractionf("voidfractionf", fvc::interpolate(voidfraction));
    surfaceScalarField rAUfvoidfraction
    (
        "(voidfraction2|A(U)F)",
        fvc::interpolate(rAUvoidfraction)
    );

    volVectorField HbyA("HbyA", rAU*UEqn.H());
    //HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            (fvc::interpolate(voidfraction*HbyA) & mesh.Sf())
          //+ rhorAUf*fvc::ddtCorr(rho, U, phi)
          + rAUfvoidfraction*fvc::ddtCorr(U, phi)
        )
    );
    fvOptions.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);

    surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
    surfaceScalarField phiGes = phiHbyA + rAUf*(fvc::interpolate(Ksl/rho)*phiS);
    if (modelType=="A")
    {
        rAUvoidfraction = volScalarField("(voidfraction2|A(U))",rAU*voidfraction*voidfraction);
    }

    //fvScalarMatrix pDDtEqn
    //(
    //    fvc::ddt(voidfraction) + voidfraction*psi*correction(fvm::ddt(p))
    //  + fvc::div(phiHbyA)
    // ==
    //    fvOptions(psi*voidfraction, p, rho.name())
    //);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            //pDDtEqn
            fvm::laplacian(rAUvoidfraction, p)
            ==
            fvc::div(phiGes)
          + particleCloud.ddtVoidfraction()
        );

        pEqn.setReference(pRefCell, pRefValue);
        //fvOptions.constrain(pEqn);

        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
        {
            pEqn.solve(mesh.solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phiGes -= pEqn.flux();
            phi = phiGes;
        }

    }
    
    //while (pimple.correctNonOrthogonal())
    //{
    //    fvScalarMatrix pEqn
    //    (
    //        //pDDtEqn
    //        fvm::laplacian(rAUvoidfraction, p)
    //        ==
    //        fvc::div(phiGes)
    //      + particleCloud.ddtVoidfraction()
    //    );

    //    //pEqn.setReference(pRefCell, pRefValue);
    //    fvOptions.constrain(pEqn);

    //    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    //    if (pimple.finalNonOrthogonalIter())
    //    {
    //        phiGes -= pEqn.flux();
    //        phi = phiGes;
    //    }

    //}
    // Second part of thermodynamic density update
    // thermo.rho() += psi*p;

    // #include "rhoEqn.H"
    // #include "compressibleContinuityErrs.H"
    #include "continuityErrorPhiPU.H"
    if (modelType == "B" || modelType == "Bfull")
        U = HbyA - rAU*fvc::grad(p) - Ksl/rho*Us*rAU;
    else
        U = HbyA - voidfraction*rAU*fvc::grad(p) - Ksl/rho*Us*rAU;

    U.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}

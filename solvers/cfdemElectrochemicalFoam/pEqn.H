{
                    volScalarField rUA = 1.0/UEqn.A();

                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

                    U = rUA*UEqn.H();

                    #ifdef version23
                    phi = ( fvc::interpolate(U*voidfraction) & mesh.Sf() )
                        + rUAfvoidfraction*fvc::ddtCorr(U, phi);
                    #else
                    phi = ( fvc::interpolate(U*voidfraction) & mesh.Sf() )
                        + fvc::ddtPhiCorr(rUAvoidfraction, U, phi);
                    #endif

                    fvOptions.makeRelative(phi);
                    adjustPhi(phi, U, p);

                    surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
                    surfaceScalarField phiGes = phi + rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                    if (modelType=="A")
                        rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                    // Non-orthogonal pressure corrector loop
                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) + particleCloud.ddtVoidfraction()
                        );
                        pEqn.setReference(pRefCell, pRefValue);

                        if
                        (
                            corr == nCorr-1
                         && nonOrth == nNonOrthCorr
                        )
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

                    } // end non-orthogonal corrector loop

                    #include "continuityErrorPhiPU.H"

                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA;
                    else
                        U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA;

    U.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}

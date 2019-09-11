{
    if(electrochemistry.active())
    {

        //volScalarField smallC
        //(
        //    IOobject
        //    (
        //        "smallC",
        //        runTime.timeName(),
        //        mesh,
        //        IOobject::NO_READ,
        //        IOobject::NO_WRITE
        //    ),
        //    mesh,
        //    dimensionedScalar("smallC", dimless, SMALL)
        //);

        // Explicitly calculate electric current density due to diffusion flux
        id *= 0.0;
        surfaceVectorField nHat = mesh.Sf()/mesh.magSf();
        forAll(Y, speciesI)
        {
            if (z[speciesI] != 0)
            {
                volVectorField gradC = fvc::grad(C[speciesI]);
                forAll(gradC.boundaryField(), patchI)
                {
                    gradC.boundaryField()[patchI] = 
                        C[speciesI].boundaryField()[patchI].snGrad()()
                        * nHat.boundaryField()[patchI];
                }
                if (transportSolver == "dilute")
                {
                    id -= F*z[speciesI]*Deff[speciesI]*gradC;
                }
                else
                {
                    id += kappa/F*t[speciesI]/z[speciesI]
                       *R*T/C[speciesI]*gradC;
                }
            }
        }
        
        // Calculate charge density
        rhoE *= 0.0;
        forAll(Y, speciesI)
        {
            if (Y[speciesI].name() != inertSpecie)
            {
                rhoE += F*z[speciesI]*C[speciesI];
            }
        }

        // Update electrochemistry model 
        electrochemistry.correct();

        label iter = 0;
        label iterMax = 1; 
        while(iter < iterMax)
        {

            // Set surface normal gradient for calculatedGradient boundary condition
            surfaceScalarField snGrad = 
                (
                    -(fvc::interpolate(il) & mesh.Sf())
                    +(fvc::interpolate(id) & mesh.Sf())
                )/(mesh.magSf()*fvc::interpolate(kappa,"interpolate(kappa)"));

            setSnGrad<calculatedGradientFvPatchScalarField>
            (
                phiEl.boundaryField(),
                snGrad.boundaryField()
            );

            // Set up matrix for Poisson equation assuming electro-neutrality
            fvScalarMatrix phiElEqn
            (
                fvm::laplacian(kappa, phiEl, "laplacian(kappa,phiE)")
                == 
               //- iTrans
                fvc::div(id) 
            );

            // Set up matrix for Poisson equation without electro-neutrality assumption
            //fvScalarMatrix phiElEqn
            //(
            //    fvm::laplacian(phiEl, "laplacian(phiE)") + rhoE/(epsilon0*epsilonR)
            //);

            // Solve for liquid electrostatic potential
            phiElEqn.solve();
            //phiEl = phiE_star*phiE_ref;

            // Explicitly calculate electric field 
            E = -fvc::grad(phiEl);
            // Explicitly calculate ionic current density
            il = kappa*E + id; 
            il.correctBoundaryConditions();
            phiEl.correctBoundaryConditions();

            iter++;
        }

        forAll(electrodes, electrodeI)
        {
            electrodes[electrodeI].correctElectricity();
        } 


    }

    // Explicitly calculate species fluxes
    forAll(Y, speciesI)
    {
        if (Y[speciesI].name() != inertSpecie
            && z[speciesI] != 0)
        {
            if (transportSolver == "dilute")
            {
                N[speciesI] = 
                   U*C[speciesI] 
                   -Deff[speciesI]*fvc::grad(C[speciesI])
                   +z[speciesI]*Deff[speciesI]/(R*T)*F*E*C[speciesI];
            }
            else
            {
                N[speciesI] = 
                   U*C[speciesI] 
                   -Deff[speciesI]*fvc::grad(C[speciesI])
                   -t[speciesI]/(z[speciesI]*F)*il;
            }
            N[speciesI].correctBoundaryConditions();
        }
    }
    forAll(electrodes, electrodeI)
    {
        electrodes[electrodeI].correctSpeciesFlux();
    } 

    if(electrochemistry.active())
    {
        // Calculate divergence of ionic current density
        div_il = fvc::div(il);
    }
}

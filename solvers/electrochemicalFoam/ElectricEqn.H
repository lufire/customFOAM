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
        forAll(Y, i)
        {
            if (z[i] != 0)
            {
                volVectorField gradC = fvc::grad(C[i]);
                forAll(gradC.boundaryField(), patchI)
                {
                    gradC.boundaryField()[patchI] = 
                        C[i].boundaryField()[patchI].snGrad()()
                        * nHat.boundaryField()[patchI];
                }
                if (transportSolver == "dilute")
                {
                    id -= F*z[i]*Deff[i]*gradC;
                }
                else
                {
                    id += kappa/F*t[i]/z[i]
                       *R*T/C[i]*gradC;
                }
            }
        }
        
        // Calculate charge density
        rhoE *= 0.0;
        forAll(Y, i)
        {
            if (Y[i].name() != inertSpecie)
            {
                rhoE += F*z[i]*C[i];
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
    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie
            && z[i] != 0)
        {
            if (transportSolver == "dilute")
            {
                N[i] = 
                   U*C[i] 
                   -Deff[i]*fvc::grad(C[i])
                   +z[i]*Deff[i]/(R*T)*F*E*C[i];
            }
            else
            {
                N[i] = 
                   U*C[i] 
                   -Deff[i]*fvc::grad(C[i])
                   -t[i]/(z[i]*F)*il;
            }
            N[i].correctBoundaryConditions();
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

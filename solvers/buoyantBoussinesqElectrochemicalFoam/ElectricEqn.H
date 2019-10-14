{
    if(electrochemistry.active())
    {
        // Explicitly calculate electric current density due to diffusion flux
        id *= 0.0;
        volScalarField divId = fvc::div(id)*0.0;
        surfaceVectorField nHat = mesh.Sf()/mesh.magSf();
        //volScalarField oneC 
        //(                
        //    IOobject
        //    (
        //        "smallC",
        //        mesh.time().timeName(),
        //        mesh,
        //        IOobject::NO_READ,
        //        IOobject::NO_WRITE
        //    ),
        //    mesh,
        //    dimensionedScalar("smallC", dimMoles/dimVol, SMALL)
        //);
        dimensionedScalar smallC("smallC", dimMoles/dimVol, SMALL);
        forAll(C, i)
        {
            if (z[i] != 0)
            {
                if (transportSolver == "dilute")
                {
                    Dm = F*z[i]*Deff[i];
                }
                else
                {
                    Dm = -kappa/F*t[i]/z[i]*R*T/(C[i]+smallC);
                }
                // Set internal field next to boundary equal to boundary values
                // (leads to better results for weighted gradient calculation)
                forAll(Dm.boundaryField(), patchI)
                {
                    fvPatchScalarField& Dmp = Dm.boundaryField()[patchI];
                    const fvPatch& patch = Dm.boundaryField()[patchI].patch();
                    const labelUList& faceCells = patch.faceCells();
                    forAll(Dmp, faceI)
                    {
                        Dm[faceCells[faceI]] = Dmp[faceI];
                    }
                }

                id -= Dm*fvc::grad(C[i], "grad(C)");
                divId -= fvc::laplacian(Dm, C[i], "laplacian(Dm,C)");
            }
        }
        
        // Calculate charge density
        rhoE *= 0.0;
        forAll(C, i)
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

        phiEl.storePrevIter();
        while(iter < iterMax)
        {
            surfaceScalarField kappas =
                fvc::interpolate(kappa, "interpolate(kappa)");
            surfaceScalarField snGrad = 
                (
                    -(fvc::interpolate(il, "interpolate(i)") & mesh.Sf())
                    +(fvc::interpolate(id, "interpolate(i)") & mesh.Sf())
                )/(mesh.magSf()*kappas);

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
                //fvc::div(id)
                divId
            );
            phiElEqn.relax();

            // Solve for liquid electrostatic potential
            phiElEqn.solve();

            phiEl.correctBoundaryConditions();
            phiEl.relax();
            // Calculate electric field 
            E = -fvc::grad(phiEl, "grad(phiE)");
           
             
            // Calculate ionic current density
            il = kappa*E + id; 
            il.correctBoundaryConditions();
            il.relax();
            //phiId -= phiElEqn.flux();

            iter++;
        }

        forAll(electrodes, electrodeI)
        {
            electrodes[electrodeI].correctElectricity();
        } 
    }

    // Explicitly calculate species fluxes
    forAll(C, i)
    {
        if (Y[i].name() != inertSpecie
            && z[i] != 0)
        {
            Di = Deff[i];
            if (transportSolver == "dilute")
            {
                N[i] = 
                   U*C[i] 
                   -Di*fvc::grad(C[i], "grad(Ci)");
                   +z[i]*mu[i]*F*E*C[i];
            }
            else
            {
                N[i] = 
                   U*C[i] 
                   -Di*fvc::grad(C[i], "grad(Ci)");
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

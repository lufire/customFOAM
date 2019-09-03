{
    if(electrochemistry.active())
    {
        volScalarField smallC
        (
            IOobject
            (
                "smallC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("smallC", dimless, SMALL)
        );

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
                //gradC.boundaryField()
                //id = F*z[speciesI]*Deff[speciesI]*fvc::grad(C[speciesI]);
                //forAll(gradC,cellI)
                //{
                //    //if(C[speciesI][cellI] > 0.0)
                //    //{
                //        id[cellI] += kappa[cellI]/F.value()*t[speciesI][cellI]/z[speciesI]
                //            *R.value()*T[cellI]/C[speciesI][cellI]*gradC[cellI];
                //    //}
                //}
                //forAll(gradC.boundaryField(),patchI)
                //{
                //    forAll(gradC.boundaryField()[patchI],faceI)
                //    {
                //        //if(C[speciesI].boundaryField()[patchI][faceI] > 0.0)
                //        //{
                //            id.boundaryField()[patchI][faceI]
                //                += kappa.boundaryField()[patchI][faceI]/F.value()
                //                    *t[speciesI].boundaryField()[patchI][faceI]
                //                    /z[speciesI]*R.value()*T.boundaryField()[patchI][faceI]
                //                    /C[speciesI].boundaryField()[patchI][faceI]
                //                    *gradC.boundaryField()[patchI][faceI];
                //        //}
                //    }
                //}
                id += kappa/F*t[speciesI]/z[speciesI]
                   *R*T/C[speciesI]*gradC;
            }
        }
        //Info << "id: " << id << endl;
        //surfaceScalarField sid = fvc::interpolate(id) & mesh.Sf();
        //tmp<volVectorField> tid = fvc::reconstruct(sid); 
        //id = tid();

        rhoE *= 0.0;
        forAll(Y, speciesI)
        {
            if (Y[speciesI].name() != inertSpecie)
            {
                rhoE += F*z[speciesI]*C[speciesI];
            }
        }
        //Info << "id interpolated: " << id << endl;
        //
        //Info << "surface diff curr density: " << sid << endl;
        //Calculate electrode charge transfer 
        electrochemistry.correct();

        //iTrans *= 0.0;
        //forAll(electrodes, electrodeI)
        //{
        //    iTrans += electrodes[electrodeI].iTrans();
        //} 
        int iter = 0;
        while(iter < 1)
        {

            // Solve for liquid electrostatic potential
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

            fvScalarMatrix phiElEqn
            (
                fvm::laplacian(kappa, phiEl, "laplacian(kappa,phiE)")
                == 
               //- iTrans
                fvc::div(id) 
            );
            //fvScalarMatrix phiElEqn
            //(
            //    fvm::laplacian(phiEl, "laplacian(phiE)") + rhoE/(epsilon0*epsilonR)
            //);

            phiElEqn.solve();
            //phiEl = phiE_star*phiE_ref;
            //const scalarField& Volume = mesh.V();
            //volScalarField phiEl1 = fvc::surfaceSum((fvc::interpolate((-il+id)/kappa) & mesh.Sf())/mesh.magSf());
            
            //phiEl.internalField() = phiEl1.internalField();//mesh.V();

            //tmp<volVectorField> tid = fvc::reconstruct(fvc::interpolate(id) & mesh.Sf()); 
            //id.internalField() = tid().internalField();
            E = -fvc::grad(phiEl);
            il = kappa*E + id; 
            ////volVectorField il2 = il*0.0; 
            il.correctBoundaryConditions();
            phiEl.correctBoundaryConditions();
            //E = (il-id)/kappa;

            iter++;
        }

        forAll(electrodes, electrodeI)
        {
            electrodes[electrodeI].correctElectricity();
        } 


        // Explicitly calculate current density and species fluxes in the electrolyte

    }

    forAll(Y, speciesI)
    {
        if (Y[speciesI].name() != inertSpecie
            && z[speciesI] != 0)
        {
            N[speciesI] = 
               U*C[speciesI] 
               -Deff[speciesI]*fvc::grad(C[speciesI])
               -t[speciesI]/(z[speciesI]*F)*il;
            N[speciesI].correctBoundaryConditions();
            //il2 += F*z[speciesI]*N[speciesI];
        }
    }
    forAll(electrodes, electrodeI)
    {
        electrodes[electrodeI].correctSpeciesFlux();
    } 
    //il = il2;
    //il.correctBoundaryConditions();

    if(electrochemistry.active())
    {
        //forAll(electrodes, electrodeI)
        //{
        //    electrodes[electrodeI].correctSpeciesFlux();
        //} 

        div_il = fvc::div(il);

    }
}

if(electrochemistry.active())
{
    // Explicitly calculate electric current density due to diffusion flux
    id *= 0.0;
    volScalarField divId = fvc::div(id)*0.0;
    surfaceVectorField nHat = mesh.Sf()/mesh.magSf();
    dimensionedScalar smallC("smallC", dimMoles/dimVol, SMALL);
    forAll(C, i)
    {
        if (z[i] != 0)
        {
            if (transportSolver == "dilute")
            {
                Dm = F*z[i]*D[i]*porosityFactor;
            }
            else
            {
                Dm = -kappa*porosityFactor*t[i]/(F*z[i])*R*T/(C[i]+smallC);
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

    //Calculate electrode charge transfer 
    electrochemistry.correct();

    volScalarField divIs = divId*0.0;
    //is *= 0.0;
    forAll(electrodes, electrodeI)
    {
        //is += electrodes[electrodeI].is();
        divIs += fvc::div(electrodes[electrodeI].is());
    }


    //iTrans *= 0.0;
    //forAll(electrodes, electrodeI)
    //{
    //    iTrans += electrodes[electrodeI].iTrans();
    //} 
    // Solve for liquid electrostatic potential
    kappaEff = kappa*porosityFactor;
    surfaceScalarField snGrad = 
        (
            -(fvc::interpolate(il) & mesh.Sf())
            +(fvc::interpolate(id) & mesh.Sf())
        )/(mesh.magSf()*fvc::interpolate(kappaEff));

    setSnGrad<calculatedGradientFvPatchScalarField>
    (
        phiEl.boundaryField(),
        snGrad.boundaryField()
    );

    Info << "Solve for electric potential in liquid electrolyte" << endl; 
    fvScalarMatrix phiElEqn
    (
        fvm::laplacian(kappaEff, phiEl, "laplacian(kappa,phiE)")
        == 
        //iTrans
        divId
      + divIs
      //  fvc::div(id) 
      //+ fvc::div(is)
    );
    phiElEqn.relax();
    phiElEqn.solve();
    phiEl.correctBoundaryConditions();

    il = kappaEff*E + id; 
    il.correctBoundaryConditions();


    forAll(electrodes, electrodeI)
    {
        electrodes[electrodeI].correctElectricity();
    } 

    // Calculate electric field
    E = -fvc::grad(phiEl);
}

// Explicitly calculate species fluxes
forAll(C, i)
{
    if (Y[i].name() != inertSpecie
        && z[i] != 0)
    {
        Di = D[i]*porosityFactor;
        if (transportSolver == "dilute")
        {
            N[i] = 
               U*C[i]*voidfraction
               -Di*fvc::grad(C[i], "grad(Ci)");
               +z[i]*mu[i]*porosityFactor*F*E*C[i];
        }
        else
        {
            N[i] = 
               U*C[i]*voidfraction 
               -Di*fvc::grad(C[i], "grad(Ci)");
               -t[i]*porosityFactor/(z[i]*F)*il;
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
 
    rhoE *= 0.0;
    forAll(Y, specieI)
    {
        if (Y[specieI].name() != inertSpecie)
        {
            rhoE += F*z[specieI]*C[specieI];
        }
    }
}

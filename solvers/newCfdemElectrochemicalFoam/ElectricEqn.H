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
    forAll(Y, specieI)
    {
        if (z[specieI] != 0)
        {
            volVectorField gradC = fvc::grad(C[specieI]);
            id = F*z[specieI]*Deff[specieI]*fvc::grad(C[specieI]);
            forAll(gradC,cellI)
            {
                if(C[specieI][cellI] > 0.0)
                {
                    id[cellI] += kappa[cellI]/F.value()*t[specieI][cellI]/z[specieI]
                        *R.value()*T[cellI]/C[specieI][cellI]*gradC[cellI];
                }
            }
            forAll(gradC.boundaryField(),patchI)
            {
                forAll(gradC.boundaryField()[patchI],faceI)
                {
                    if(C[specieI].boundaryField()[patchI][faceI] > 0.0)
                    {
                        id.boundaryField()[patchI][faceI]
                            += kappa.boundaryField()[patchI][faceI]/F.value()
                                *t[specieI].boundaryField()[patchI][faceI]
                                /z[specieI]*R.value()*T.boundaryField()[patchI][faceI]
                                /C[specieI].boundaryField()[patchI][faceI]
                                *gradC.boundaryField()[patchI][faceI];
                    }
                }
            }
            //id += kappa/F*t[specieI]/z[specieI]
            //   *R*T/(C[specieI]+smallC)*fvc::grad(C[specieI]);
        }
    }

    //Calculate electrode charge transfer 
    electrochemistry.correct();

    is *= 0.0;
    forAll(electrodes, electrodeI)
    {
        is += electrodes[electrodeI].is();
    }

    // Calculate conductivity field
    //kappa *= 0;         
    //forAll(C, specieI)
    //{
    //    //C[specieI].correctBoundaryConditions();
    //    if (Y[specieI].name() != inertSpecie)
    //    {
    //        kappa += F*F*z[specieI]*z[specieI]*mu[specieI]*C[specieI];
    //    }
    //}

    //// Explicitly calculate electric current density due to diffusion flux
    //id *= 0.0;
    //forAll(Y, specieI)
    //{
    //    if (Y[specieI].name() != inertSpecie)
    //    {
    //        id -= F*z[specieI]*DEff[specieI]*fvc::grad(C[specieI]);
    //    }
    //}
    // Correct diffusion current density due to voidfraction (Bruggeman)
    //id *= pow(voidfraction,1.5);

    //surfaceScalarField idf(fvc::interpolate(id)&mesh.Sf());
    //rhoE *= 0.0;
    kappaeff = kappa*pow(voidfraction,1.5);

    iTrans *= 0.0;
    forAll(electrodes, electrodeI)
    {
        iTrans += electrodes[electrodeI].iTrans();
    } 
    // Solve for liquid electrostatic potential
    surfaceScalarField snGrad = 
        (
            -(fvc::interpolate(il) & mesh.Sf())
            +(fvc::interpolate(id) & mesh.Sf())
        )/(mesh.magSf()*fvc::interpolate(kappaeff,"interpolate(kappa)"));

    setSnGrad<calculatedGradientFvPatchScalarField>
    (
        phiEl.boundaryField(),
        snGrad.boundaryField()
    );
    Info << "Solve for electric potential in liquid electrolyte" << endl; 
    fvScalarMatrix phiElEqn
    (
        fvm::laplacian(kappaeff, phiEl, "laplacian(kappa,phiEl)")
      //+ fvc::div(id)
      //+ fvc::div(is)
        == 
        //iTrans
        fvc::div(id) 
      + fvc::div(is)
    );
    phiElEqn.solve();
    phiEl.correctBoundaryConditions();

    il = kappaeff*E + id; 
    //volVectorField il2 = il*0.0; 
    il.correctBoundaryConditions();
    // Solve Poisson equation
    //dimensionedScalar epsilon = epsilon0*epsilonR;
    //fvScalarMatrix phiElEqn
    //(
    //    fvm::laplacian(epsilon, phiEl, "laplacian(epsilon, phiEl)") + rhoE
    //);


    forAll(electrodes, electrodeI)
    {
        electrodes[electrodeI].correctElectricity();
    } 
    //phiEl = phiEl*(1.0-0.7)+phiEl.oldTime()*(0.7);
    //phiEl *= 0;

    // Calculate electric field
    E = -fvc::grad(phiEl);
}

forAll(Y, specieI)
{
    if (Y[specieI].name() != inertSpecie
        && z[specieI] != 0)
    {
        N[specieI] = 
           U*C[specieI]*voidfraction
           -Deff[specieI]*fvc::grad(C[specieI])
           -teff[specieI]/(z[specieI]*F)*il;
        N[specieI].correctBoundaryConditions();
        //il2 += F*z[specieI]*N[specieI];
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
 
    // Calculate liquid current density
    //il = id + kappa*E;
    //im = kappa*E;
    //il.correctBoundaryConditions();


    rhoE *= 0.0;
    forAll(Y, specieI)
    {
        if (Y[specieI].name() != inertSpecie)
        {
            rhoE += F*z[specieI]*C[specieI];
        }
    }
    //div_il = fvc::div(il);
    //div_is = fvc::div(is);
}
}

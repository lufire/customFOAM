{
    // Update effective diffusion coefficients 
    forAll(D, i)
    {
        Deff[i] = D[i];
    }
    //reaction->correct();
    //dQ = reaction->dQ();
    label inertIndex = -1;
    label balanceIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    zC *= 0.0;
    //il *= 0.0;
    
    //volScalarField iSource = electrodes[0].stoichCoeff(0)
    //    /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();

    forAll(C, i)
    {

        //volScalarField& Yi = Y[i];       
        //if (Y[i].name() != inertSpecie)
        if (Y[i].name() == electrolyte.inertSpecies())
        {
            inertIndex = i;
        }
        else if (Y[i].name() == electrolyte.balanceSpecies())
        {
            balanceIndex = i;
        }
        else
        {
            volScalarField& Ci = C[i];
            if (transportSolver == "dilute")
            {
                mu[i] = Deff[i]/(R*T);
                surfaceScalarField zmuFEf 
                (
                    "zmuFEf", 
                    (
                        (fvc::interpolate(z[i]*mu[i]*F*E, "interpolate(zmuFE)")
                         & mesh.Sf())
                    )
                );

                // Set surface normal gradient for boundary field
                surfaceScalarField snGrad = 
                    (
                        phi*fvc::interpolate(Ci)
                        -(mesh.Sf() & fvc::interpolate(N[i]))
                        -(zmuFEf * fvc::interpolate(Ci))
                    )/(mesh.magSf()*fvc::interpolate(Deff[i], "interpolate(D)"));
                
                setSnGrad<calculatedGradientFvPatchScalarField>
                (
                    Ci.boundaryField(),
                    snGrad.boundaryField()
                );

                // Species transport equation
                fvScalarMatrix CiEqn
                (
                    fvm::ddt(Ci)
                  + fvm::div(phi, Ci) 
                  - fvm::laplacian(Deff[i],Ci, "laplacian(D,Ci)")
                  + fvm::div(zmuFEf,Ci, "div(zmuFEf,Ci)")
                    ==
                    fvOptions(Ci)
                );
                CiEqn.relax();
                fvOptions.constrain(CiEqn);
                CiEqn.solve(mesh.solver("Ci"));
            }
            else
            {
                // Set surface normal gradient for boundary field
                surfaceScalarField snGrad = 
                    (
                        phi*fvc::interpolate(Ci)
                        -(mesh.Sf() & fvc::interpolate(N[i]))
                        -(mesh.Sf() 
                         & fvc::interpolate(t[i]/(z[i]*F), "interpolate(t)")
                             *fvc::interpolate(il))
                    )/(mesh.magSf()*fvc::interpolate(Deff[i], "interpolate(D)"));
                
                setSnGrad<calculatedGradientFvPatchScalarField>
                (
                    Ci.boundaryField(),
                    snGrad.boundaryField()
                );

                // Species transport equation
                fvScalarMatrix CiEqn
                (
                    fvm::ddt(Ci)
                  + fvm::div(phi, Ci) 
                  - fvm::laplacian(Deff[i],Ci, "laplacian(D,Ci)")
                    ==
                    fvc::div(t[i]/(z[i]*F)*il)
                  + fvOptions(Ci)
                );
                CiEqn.relax();
                fvOptions.constrain(CiEqn);
                CiEqn.solve(mesh.solver("Ci"));
            }


            Ci.correctBoundaryConditions();
            fvOptions.correct(Ci);
            Ci.max(SMALL);

            zC += z[i]*Ci;
            
            // Calculate species mass fraction
            Y[i] = Ci*M[i]/rho;
            Yt += Y[i];
        }
        //if (Y[i].name() != inertSpecie
        //    && z[i] != 0
        //    && chargedSpeciesCount == nChargedSpecies)
        //{
        //    volScalarField& Ci = C[i];
        //    Ci = -zC/z[i];
        //    Ci.correctBoundaryConditions();
        //    fvOptions.correct(Ci);
        //    Ci.max(0.0);

        //    // Calculate species mass fraction
        //    Y[i] = Ci*M[i]/rho;
        //    Yt += Y[i];
        //}
    }

    forAll(electrodes, electrodeI)
    {
        electrodes[electrodeI].correctSpecies(C);
    }       

    { 
        volScalarField& Ci = C[balanceIndex];
        Ci = -zC/z[balanceIndex];
        Ci.correctBoundaryConditions();
        fvOptions.correct(Ci);
        Ci.max(0.0);
        // Calculate species mass fraction
        Y[balanceIndex] = Ci*M[balanceIndex]/rho;
        Yt += Y[balanceIndex];
    }
  
    //forAll(C,i)
    //{
    //    if (Y[i].name() != inertSpecie
    //        && z[i] != 0)
    //    {
    //        N[i] = 
    //           U*C[i] 
    //           -Deff[i]*fvc::grad(C[i])
    //           -t[i]/(z[i]*F)*il;
    //        N[i].correctBoundaryConditions();
    //    }
    //}
    //forAll(electrodes, electrodeI)
    //{
    //    electrodes[electrodeI].correctSpeciesFlux();
    //} 
    Y[inertIndex] = scalar(1) - Yt;
    C[inertIndex] = Y[inertIndex]*rho/M[inertIndex];
    C[inertIndex].max(0.0);
    Y[inertIndex].max(0.0);

    //mstm().updateMolarFractions();
}

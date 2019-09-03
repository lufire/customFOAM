{
    // Update effective diffusion coefficients 
    forAll(D, speciesI)
    {
        Deff[speciesI] = D[speciesI];
        //mu[speciesI] = Deff[speciesI]/(R*T);
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

    forAll(C, speciesI)
    {

        //volScalarField& Yi = Y[speciesI];       
        //if (Y[speciesI].name() != inertSpecie)
        if (Y[speciesI].name() == electrolyte.inertSpecies())
        {
            inertIndex = speciesI;
        }
        else if (Y[speciesI].name() == electrolyte.balanceSpecies())
        {
            balanceIndex = speciesI;
        }
        else
        {

            //dimensionedScalar Pe = u_ref*L_ref/Deff[speciesI];
            
            surfaceScalarField snGrad = 
                (
                    phi*fvc::interpolate(C[speciesI])
                    -(mesh.Sf() & fvc::interpolate(N[speciesI]))
                    -(mesh.Sf() 
                     & fvc::interpolate(
                            t[speciesI]/(z[speciesI]*F),
                            "interpolate(t)")
                         *fvc::interpolate(il))
                )/(mesh.magSf()*fvc::interpolate(
                        Deff[speciesI],
                        "interpolate(D)"));
            
            setSnGrad<calculatedGradientFvPatchScalarField>
            (
                C[speciesI].boundaryField(),
                snGrad.boundaryField()
            );

            volScalarField& Ci = C[speciesI];
            fvScalarMatrix CiEqn
            (
                fvm::ddt(Ci)
              + fvm::div(phi, Ci) 
              //+ mvConvection->fvmDiv(phi, Ci) 
              //- fvm::Sp(fvc::div(phi), Ci)
              - fvm::laplacian
                (
                    Deff[speciesI],Ci, "laplacian(D,Ci)"
                )
              //+ fvm::div(zmuFEf,Ci, "div(zmuFEf,Ci)")
                ==
              - fvc::div(t[speciesI]/(z[speciesI]*F)*il)
              //- fvc::div(Ni)
              //- iSource
              + fvOptions(Ci)
              //  gamma/M[speciesI]
              //  *mstm().RR(Y[speciesI].name(), chemistry, speciesI)
            );
            //Info << "iSource[0]: " << iSource[0] << endl;
            //Info << "iSource[1]: " << iSource[1] << endl;
            CiEqn.relax();

            fvOptions.constrain(CiEqn);

            CiEqn.solve(mesh.solver("Ci"));
            //C[speciesI] = Ci*C_ref;

            Ci.correctBoundaryConditions();


            fvOptions.correct(Ci);
            Ci.max(SMALL);
            

            zC += z[speciesI]*Ci;
            
            // Calculate species mass fraction
            Y[speciesI] = Ci*M[speciesI]/rho;
            Yt += Y[speciesI];
        }
        //if (Y[speciesI].name() != inertSpecie
        //    && z[speciesI] != 0
        //    && chargedSpeciesCount == nChargedSpecies)
        //{
        //    volScalarField& Ci = C[speciesI];
        //    Ci = -zC/z[speciesI];
        //    Ci.correctBoundaryConditions();
        //    fvOptions.correct(Ci);
        //    Ci.max(0.0);

        //    // Calculate species mass fraction
        //    Y[speciesI] = Ci*M[speciesI]/rho;
        //    Yt += Y[speciesI];
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
  
    //forAll(C,speciesI)
    //{
    //    if (Y[speciesI].name() != inertSpecie
    //        && z[speciesI] != 0)
    //    {
    //        N[speciesI] = 
    //           U*C[speciesI] 
    //           -Deff[speciesI]*fvc::grad(C[speciesI])
    //           -t[speciesI]/(z[speciesI]*F)*il;
    //        N[speciesI].correctBoundaryConditions();
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

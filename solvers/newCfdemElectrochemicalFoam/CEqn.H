{
    //mstm().updateCoefficients();
   
    //Update effective diffusion coefficients 
    forAll(D, specieI)
    {
        Deff[specieI] = D[specieI]*pow(voidfraction,1.5);
        //mu[specieI] = Deff[specieI]*pow(voidfraction,1.5);
        teff[specieI] = t[specieI]*pow(voidfraction,1.5);
    }
    //reaction->correct();
    //dQ = reaction->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    zC *= 0.0;
    //il *= 0.0;
    
    label chargedSpeciesCount = 0;

    volScalarField iSource = electrodes[0].stoichCoeff(0)
        /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();

    forAll(C, specieI)
    {

        volScalarField& Ci = C[specieI];
        volVectorField& Ni = N[specieI];
        //volScalarField& Yi = Y[specieI];       
        if (z[specieI] != 0)
        {
            chargedSpeciesCount++;
        }
//        if (Y[specieI].name() != inertSpecie)
        if (Y[specieI].name() != inertSpecie
            && chargedSpeciesCount < nChargedSpecies)
        {
            
            //surfaceScalarField zmuFEf 
            //(
            //    "zmuFEf", 
            //    (
            //        (fvc::interpolate(z[specieI]*mu[specieI]*F*E) 
            //         & mesh.Sf())
            //    )//*fvc::interpolate(z(specieI)*mu(spec*Di/(T*kB))
            //);
            
            iSource = electrodes[0].stoichCoeff(specieI)
                /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();
            for(label i=1; i<electrodes.size(); ++i)
            {
                iSource += electrodes[i].stoichCoeff(specieI)
                    /(electrodes[i].electronNumber()*F)*electrodes[i].iTrans();
            }
            //Info << " Name of Species: " << Ci.name() << endl;
            //Info << " iSource[0]: " << iSource[0] << endl;
            surfaceScalarField snGrad = 
                (
                    phi*fvc::interpolate(voidfraction*Ci)
                    -(mesh.Sf() & fvc::interpolate(Ni))
                    -(mesh.Sf() 
                     & fvc::interpolate(
                            teff[specieI]/(z[specieI]*F),
                            "interpolate(t)")
                         *fvc::interpolate(il))
                )/(mesh.magSf()*fvc::interpolate(
                        Deff[specieI],
                        "interpolate(D)"));
            
            setSnGrad<calculatedGradientFvPatchScalarField>
            (
                Ci.boundaryField(),
                snGrad.boundaryField()
            );
            
            fvScalarMatrix CiEqn
            (
                fvm::ddt(voidfraction,Ci) - fvm::Sp(fvc::ddt(voidfraction),Ci)
              + fvm::div(phi, Ci) - fvm::Sp(fvc::div(phi), Ci)
              //+ mvConvection->fvmDiv(phi, Ci) 
              //- fvm::Sp(fvc::div(phi), Ci)
              - fvm::laplacian
                (
                    Deff[specieI],Ci, "laplacian(D,Ci)"
                )
              //+ fvm::div(zmuFEf,Ci, "div(zmuFEf,Ci)")
                ==
              - fvc::div(teff[specieI]/(z[specieI]*F)*il)
              - iSource
              + fvOptions(Ci)
              //  gamma/M[specieI]
              //  *mstm().RR(Y[specieI].name(), chemistry, specieI)
            );

            CiEqn.relax();

            fvOptions.constrain(CiEqn);

            CiEqn.solve(mesh.solver("Ci"));

            fvOptions.correct(Ci);
            
            Ci.max(0.0);

            //Ni = -Deff[specieI]*fvc::grad(Ci) + z[specieI]*mu[specieI]*F*E*Ci + U*Ci;
            //Ni.correctBoundaryConditions();

            Ci.correctBoundaryConditions();

            zC += z[specieI]*Ci;

            //il += F*z[specieI]*Ni;
            
            // Calculate species mass fraction
            Y[specieI] = Ci*M[specieI]/rho;
            Yt += Y[specieI];
        }
        else if (Y[specieI].name() != inertSpecie
            && z[specieI] != 0
            && chargedSpeciesCount == nChargedSpecies)
        {
            Ci = -zC/z[specieI];

            Ci.correctBoundaryConditions();
            fvOptions.correct(Ci);
            Ci.max(0.0);

            //Ni = -Deff[specieI]*fvc::grad(Ci) + z[specieI]*mu[specieI]*F*E*Ci;
            // Calculate species flux density
            //N[specieI] = z[specieI]*mu[specieI]*F*Ci*E-D[specieI]*fvc::grad(Ci);

            //il += F*z[specieI]*Ni;
            // Calculate species mass fraction
            Y[specieI] = Ci*M[specieI]/rho;
            Yt += Y[specieI];
        }
        else if (Y[specieI].name() == inertSpecie)
        {
            inertIndex = specieI;
        }
    }

    forAll(electrodes, electrodeI)
    {
        electrodes[electrodeI].correctSpecies(C);
    } 

    Y[inertIndex] = scalar(1) - Yt;
    C[inertIndex] = Y[inertIndex]*rho/M[inertIndex];
    C[inertIndex].max(0.0);
    Y[inertIndex].max(0.0);
    //il.correctBoundaryConditions();

    //mstm().updateMolarFractions();
}

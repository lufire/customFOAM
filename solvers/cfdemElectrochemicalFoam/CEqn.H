{
    //reaction->correct();
    //dQ = reaction->dQ();
    label balanceIndex = -1;
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    zC *= 0.0;

    //Update effective transport properties 
    forAll(D, i)
    {
        //Deff[i] = D[i]*pow(voidfraction,1.5);
        //tEff[i] = t[i]*pow(voidfraction,1.5);
        mu[i] = D[i]/(R*T);
    }
    //label chargedSpeciesCount = 0;

    volScalarField iSource = electrodes[0].stoichCoeff(0)
        /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();

    forAll(C, i)
    {

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
            volVectorField& Ni = N[i];

            volScalarField Di = D[i]*porosityFactor;
            
            // Sum up source terms due to electrochemical reactions
            iSource = electrodes[0].stoichCoeff(i)
                /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();
            for(label j=1; j<electrodes.size(); ++j)
            {
                iSource += electrodes[j].stoichCoeff(i)
                    /(electrodes[j].electronNumber()*F)*electrodes[j].iTrans();
            }

            if (transportSolver == "dilute")
            {
                surfaceScalarField zmuFEf 
                (
                    "zmuFEf", 
                    (
                        (fvc::interpolate(porosityFactor*z[i]*mu[i]*F*E) 
                         & mesh.Sf())
                    )
                );
                surfaceScalarField snGrad = 
                    (
                        phi*fvc::interpolate(voidfraction*Ci)
                        -(mesh.Sf() & fvc::interpolate(Ni))
                        +(zmuFEf * fvc::interpolate(Ci))
                    )/(mesh.magSf()*fvc::interpolate(Di));
                fvScalarMatrix CiEqn
                (
                    fvm::ddt(voidfraction,Ci) - fvm::Sp(fvc::ddt(voidfraction),Ci)
                  + fvm::div(phi, Ci) - fvm::Sp(fvc::div(phi), Ci)
                  - fvm::laplacian(Di, Ci, "laplacian(D,C)")
                  //+ fvm::div(zmuFEf,Ci, "div(zmuFEf,C)")
                  - fvc::laplacian(porosityFactor*z[i]*mu[i]*F*Ci, phiEl, 
                                   "laplacian(zmuFC,phiE)")
                    ==
                  - iSource
                  + fvOptions(Ci)
                  //  gamma/M[i]
                  //  *mstm().RR(Y[i].name(), chemistry, i)
                );
                CiEqn.relax();
                fvOptions.constrain(CiEqn);
                CiEqn.solve(mesh.solver("C"));
            }
            else
            {
                volScalarField ti = porosityFactor*t[i];
                surfaceScalarField snGrad = 
                    (
                        phi*fvc::interpolate(voidfraction*Ci)
                        -(mesh.Sf() & fvc::interpolate(Ni))
                        -(mesh.Sf() & fvc::interpolate(ti*il/(z[i]*F)))
                    )/(mesh.magSf()*fvc::interpolate(Di));
                
                setSnGrad<calculatedGradientFvPatchScalarField>
                (
                    Ci.boundaryField(),
                    snGrad.boundaryField()
                );
                
                fvScalarMatrix CiEqn
                (
                    fvm::ddt(voidfraction,Ci) - fvm::Sp(fvc::ddt(voidfraction),Ci)
                  + fvm::div(phi, Ci) - fvm::Sp(fvc::div(phi), Ci)
                  - fvm::laplacian(Di, Ci, "laplacian(D,C)")
                    ==
                  - fvc::div(ti/(z[i]*F)*il)
                  - iSource
                  + fvOptions(Ci)
                  //  gamma/M[i]
                  //  *mstm().RR(Y[i].name(), chemistry, i)
                );

                CiEqn.relax();
                fvOptions.constrain(CiEqn);
                CiEqn.solve(mesh.solver("C"));
            }
            fvOptions.correct(Ci);
            Ci.max(0.0);

            //Ni = -Deff[i]*fvc::grad(Ci) + z[i]*mu[i]*F*E*Ci + U*Ci;
            //Ni.correctBoundaryConditions();

            Ci.correctBoundaryConditions();

            zC += z[i]*Ci;

            //il += F*z[i]*Ni;
            
            // Calculate species mass fraction
            Y[i] = Ci*M[i]/rho;
            Yt += Y[i];
        }
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
        Ci.max(SMALL);
        // Calculate species mass fraction
        Y[balanceIndex] = Ci*M[balanceIndex]/rho;
        Yt += Y[balanceIndex];
    }

    Y[inertIndex] = scalar(1) - Yt;
    C[inertIndex] = Y[inertIndex]*rho/M[inertIndex];
    C[inertIndex].max(0.0);
    Y[inertIndex].max(0.0);
    //il.correctBoundaryConditions();

    //mstm().updateMolarFractions();
}

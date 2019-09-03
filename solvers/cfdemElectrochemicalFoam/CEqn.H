{
    mstm().updateCoefficients();
   
    //Update effective diffusion coefficients 
    forAll(D, speciesI)
    {
        DEff[speciesI] = D[speciesI]*pow(voidfraction,1.5);
        mu[speciesI] = DEff[speciesI]*pow(voidfraction,1.5);
    }
    //reaction->correct();
    //dQ = reaction->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    zC *= 0.0;
    il *= 0.0;
    
    label chargedSpeciesCount = 0;

    volScalarField iSource = electrodes[0].stoichCoeff(0)
        /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();

    forAll(C, speciesI)
    {

        volScalarField& Ci = C[speciesI];
        volVectorField& Ni = N[speciesI];
        //volScalarField& Yi = Y[speciesI];       
        if (z[speciesI] != 0)
        {
            chargedSpeciesCount++;
        }
//        if (Y[speciesI].name() != inertSpecie)
        if (Y[speciesI].name() != inertSpecie
            && chargedSpeciesCount < nChargedSpecies)
        {
            
            surfaceScalarField zmuFEf 
            (
                "zmuFEf", 
                (
                    (fvc::interpolate(z[speciesI]*mu[speciesI]*F*E) 
                     & mesh.Sf())
                )//*fvc::interpolate(z(speciesI)*mu(spec*Di/(T*kB))
            );
            
            iSource = electrodes[0].stoichCoeff(speciesI)
                /(electrodes[0].electronNumber()*F)*electrodes[0].iTrans();
            for(label i=1; i<electrodes.size(); ++i)
            {
                iSource += electrodes[i].stoichCoeff(speciesI)
                    /(electrodes[i].electronNumber()*F)*electrodes[i].iTrans();
            }
            //Info << " Name of Species: " << Ci.name() << endl;
            //Info << " iSource[0]: " << iSource[0] << endl;
            
            fvScalarMatrix CiEqn
            (
                fvm::ddt(voidfraction,Ci) - fvm::Sp(fvc::ddt(voidfraction),Ci)
              + fvm::div(phi, Ci) - fvm::Sp(fvc::div(phi), Ci)
              //+ mvConvection->fvmDiv(phi, Ci) 
              //- fvm::Sp(fvc::div(phi), Ci)
              - fvm::laplacian
                (
                    DEff[speciesI],Ci, "laplacian(D,Ci)"
                )
              + fvm::div(zmuFEf,Ci, "div(zmuFEf,Ci)")
                ==
              - iSource
              + fvOptions(Ci)
              //  gamma/M[speciesI]
              //  *mstm().RR(Y[speciesI].name(), chemistry, speciesI)
            );

            CiEqn.relax();

            fvOptions.constrain(CiEqn);

            CiEqn.solve(mesh.solver("Ci"));

            fvOptions.correct(Ci);
            
            Ci.max(0.0);

            Ni = -DEff[speciesI]*fvc::grad(Ci) + z[speciesI]*mu[speciesI]*F*E*Ci + U*Ci;
            Ni.correctBoundaryConditions();

            Ci.correctBoundaryConditions();

            zC += z[speciesI]*Ci;
            //forAll(zC.boundaryField(), patchI)
            //{
            //    forAll(zC.boundaryField()[patchI], faceI)
            //    {
            //        zC.boundaryField()[patchI][faceI] +=  Ci.boundaryField()[patchI][faceI]*z[speciesI];
            //        if (zC.boundaryField()[patchI].patch().name()== "currentCollector")
            //        {
            //            Info << " Name of Species: " << Ci.name() << endl;
            //            Info << " Ci = " << Ci.boundaryField()[patchI][faceI] << "; zC = " << zC.boundaryField()[patchI][faceI] << endl;
            //        }
            //        if (zC.boundaryField()[patchI].patch().name()== "separator")
            //        {
            //            Info << " Name of Species: " << Ci.name() << endl;
            //            Info << " Ci = " << Ci.boundaryField()[patchI][faceI] << "; zC = " << zC.boundaryField()[patchI][faceI] << endl;
            //        }
            //    }
            //}
            il += F*z[speciesI]*Ni;
            
            // Calculate species mass fraction
            Y[speciesI] = Ci*M[speciesI]/rho;
            Yt += Y[speciesI];
        }
        else if (Y[speciesI].name() != inertSpecie
            && z[speciesI] != 0
            && chargedSpeciesCount == nChargedSpecies)
        {
            Ci = -zC/z[speciesI];
            //forAll(Ci.boundaryField(), patchI)
            //{
            //    forAll(Ci.boundaryField()[patchI], faceI)
            //    {
            //        //Ci.boundaryField()[patchI][faceI] =  -zC.boundaryField()[patchI][faceI]/z[speciesI];
            //        if (zC.boundaryField()[patchI].patch().name()== "currentCollector")
            //        {
            //            Info << " Name of Species: " << Ci.name() << endl;
            //            Info << " Ci = " << Ci.boundaryField()[patchI][faceI] << "; zC = " << zC.boundaryField()[patchI][faceI] << endl;
            //        }
            //        if (zC.boundaryField()[patchI].patch().name()== "separator")
            //        {
            //            Info << " Name of Species: " << Ci.name() << endl;
            //            Info << " Ci = " << Ci.boundaryField()[patchI][faceI] << "; zC = " << zC.boundaryField()[patchI][faceI] << endl;
            //        }
            //    }
            //}
            Ci.correctBoundaryConditions();
            fvOptions.correct(Ci);
            Ci.max(0.0);

            Ni = -DEff[speciesI]*fvc::grad(Ci) + z[speciesI]*mu[speciesI]*F*E*Ci;
            // Calculate species flux density
            //N[speciesI] = z[speciesI]*mu[speciesI]*F*Ci*E-D[speciesI]*fvc::grad(Ci);

            il += F*z[speciesI]*Ni;
            // Calculate species mass fraction
            Y[speciesI] = Ci*M[speciesI]/rho;
            Yt += Y[speciesI];
        }
        else if (Y[speciesI].name() == inertSpecie)
        {
            inertIndex = speciesI;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    C[inertIndex] = Y[inertIndex]*rho/M[inertIndex];
    C[inertIndex].max(0.0);
    Y[inertIndex].max(0.0);
    il.correctBoundaryConditions();

    mstm().updateMolarFractions();
}

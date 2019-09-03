tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


// A.Alexiou 2015
// mstm().correct(kappa, chemistry, fields); // OF 2.1
// ------------------------------------------------------------------------- //
// A.Alexiou 2015
//mstm().correct(Ve, rhoE, Y, kappa, chemistry, fields); // OF 2.3
// ------------------------------------------------------------------------- //

{
    reaction->correct();
    mstm().updateCoefficients();
    dQ = reaction->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    volScalarField zC(Yt*mstm().z(0)/mstm().W(0)
            *dimensionedScalar("1000",dimensionSet(-2, 0, 0, 0, 1, 0, 0),1000));
    
    label chargedSpeciesCount = 0;
    forAll(Y, speciesI)
    {
        volScalarField& Yi = Y[speciesI];    
        if (mstm().z(speciesI) != 0)
        {
            chargedSpeciesCount++;
        }
        if (Y[speciesI].name() != inertSpecie)
//        if (Y[speciesI].name() != inertSpecie 
//            && chargedSpeciesCount < nChargedSpecies)
        {


            const volScalarField& Di = mstm().Di(speciesI);
            
            surfaceScalarField zeMEf 
            (
                "zeMEf", 
                (
                    (fvc::interpolate(E) & mesh.Sf())              
                )*fvc::interpolate(mstm().z(speciesI)*e*Di/(T*kB))
            );                

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
//              - fvm::laplacian(turbulence->muEff(), Yi)
              - fvm::laplacian(Di,Yi, "laplacian(D,Yi)")
              + fvm::div(zeMEf,Yi, "div(zeMEf,Yi)")
             ==
                kappa*mstm().RR(Yi.name(), chemistry, speciesI)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));
            //YiEqn.solve();

            fvOptions.correct(Yi);

            Yi.max(0.0);
            zC += Yi*mstm().z(speciesI)/mstm().W(speciesI)
               *dimensionedScalar(
                    "1000",dimensionSet(-2, 0, 0, 0, 1, 0, 0),1000);
            Yt += Yi;

            if(Y[speciesI].name() == anodeReductantName)
            {
                anodeReductant = Y[speciesI];
            }
            else if(Y[speciesI].name() == anodeOxidantName)
            {
                anodeOxidant = Y[speciesI];
            }
            else if(Y[speciesI].name() == cathodeReductantName)
            {
                cathodeReductant = Y[speciesI];
            }
            else if(Y[speciesI].name() == cathodeOxidantName)
            {
                cathodeOxidant = Y[speciesI];
            }                        
        }
/*        else if (Y[speciesI].name() != inertSpecie
            && mstm().z(speciesI) != 0
            && chargedSpeciesCount == nChargedSpecies)
        {
            Yi = -mstm().W(speciesI)/dimensionedScalar(
                "1000",dimensionSet(-2, 0, 0, 0, 1, 0, 0),1000)
               /mstm().z(speciesI)*zC;
            //Yi.correctBoundaryConditions();
            //fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }*/
        else if (Y[speciesI].name() == inertSpecie)
        {
            inertIndex = speciesI;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    mstm().updateMolarFractions();
}

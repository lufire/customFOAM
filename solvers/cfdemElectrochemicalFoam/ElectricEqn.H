//Calculate electrode charge transfer 
electrochemistry.correct();

is *= 0.0;
forAll(electrodes, electrodeI)
{
    is += electrodes[electrodeI].is();
}
//Update concentration, diffusion coefficients and mobility
//forAll(Y, speciesI)
//{
//    //C[speciesI] = Y[speciesI]*rho/M[speciesI];
//    if (Y[speciesI].name() != inertSpecie)
//    {    
//        mu[speciesI] = DEff[speciesI]/(R*T);
//    }
//}

// Calculate conductivity field
kappa *= 0;         
forAll(C, speciesI)
{
    //C[speciesI].correctBoundaryConditions();
    if (Y[speciesI].name() != inertSpecie)
    {
        kappa += F*F*z[speciesI]*z[speciesI]*mu[speciesI]*C[speciesI];
    }
}
//surfaceScalarField kappaf("kappaf",fvc::interpolate(kappa));
//dimensionedScalar kappaInit
//(
//    "kappa", 
//    sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
//    1.0e4 
//);
//kappa += kappaInit*Foam::exp(-100.0*runTime.timeOutputValue());


// Explicitly calculate electric current density due to diffusion flux
id *= 0.0;
forAll(Y, speciesI)
{
    if (Y[speciesI].name() != inertSpecie)
    {
        id -= F*z[speciesI]*DEff[speciesI]*fvc::grad(C[speciesI]);
    }
}
// Correct diffusion current density due to voidfraction (Bruggeman)
//id *= pow(voidfraction,1.5);

//surfaceScalarField idf(fvc::interpolate(id)&mesh.Sf());
rhoE *= 0.0;
forAll(Y, speciesI)
{
    if (Y[speciesI].name() != inertSpecie)
    {
        rhoE += F*z[speciesI]*C[speciesI];
    }
}
//rhoE *= 0.0;

iTrans *= 0.0;
forAll(electrodes, electrodeI)
{
    iTrans += electrodes[electrodeI].iTrans();
} 
// Solve for liquid electrostatic potential
Info << "Solve for electric potential in liquid electrolyte" << endl; 
fvScalarMatrix phiElEqn
(
    fvm::laplacian(kappa, phiEl, "laplacian(kappa,phiEl)")
  //+ fvc::div(id)
  //+ fvc::div(is)
    == 
    //iTrans
    fvc::div(id) 
  + fvc::div(is)
);


// Solve Poisson equation
//dimensionedScalar epsilon = epsilon0*epsilonR;
//fvScalarMatrix phiElEqn
//(
//    fvm::laplacian(epsilon, phiEl, "laplacian(epsilon, phiEl)") + rhoE
//);

phiElEqn.solve();
phiEl.correctBoundaryConditions();

//phiEl = phiEl*(1.0-0.7)+phiEl.oldTime()*(0.7);
//phiEl *= 0;

// Calculate electric field
E = -fvc::grad(phiEl);

// Calculate liquid current density
//il = id + kappa*E;
//im = kappa*E;
//il.correctBoundaryConditions();



//div_il = fvc::div(il);
//div_is = fvc::div(is);


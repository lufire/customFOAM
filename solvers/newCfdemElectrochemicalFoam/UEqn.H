    fvVectorMatrix UEqn
    (
        fvm::ddt(voidfraction, U) - fvm::Sp(fvc::ddt(voidfraction),U)
      + fvm::div(phi, U) - fvm::Sp(fvc::div(phi),U)
      //+ turbulence->divDevRhoReff(U)
      + particleCloud.divVoidfractionTau(U, voidfraction)
     ==
      - fvm::Sp(Ksl/rho,U)
      //  rho*g
      + fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    //if (pimple.momentumPredictor())
    //{
    if (momentumPredictor && (modelType=="B" || modelType=="Bfull"))
    {
        solve(UEqn == - fvc::grad(p) + Ksl/rho*Us);
        fvOptions.correct(U);
    }
    else if (momentumPredictor)
    {
        solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us);
        fvOptions.correct(U);
    }

    K = 0.5*magSqr(U);

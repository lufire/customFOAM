// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    //fvm::ddt(U)
    fvm::div(phi, U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

//volScalarField rAU(1.0/UEqn().A());

if (simple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}

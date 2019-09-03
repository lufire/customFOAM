// Solve the Momentum equation
tmp<fvVectorMatrix> UEqn
(
    fvm::div(phi, U)
  //- fvm::laplacian(nu, U)
  //- fvc::div(nu*dev(T(fvc::grad(U))), "div(nu*dev(T(grad(U))))") 
  - rheology.divTau(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if (simple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
}

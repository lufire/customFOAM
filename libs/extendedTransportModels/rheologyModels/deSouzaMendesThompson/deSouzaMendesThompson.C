/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more dnuils.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "deSouzaMendesThompson.H"
#include "addToRunTimeSelectionTable.H"
#include "solutionControl.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(deSouzaMendesThompson, 0);
addToRunTimeSelectionTable(rheologyModel, deSouzaMendesThompson, dictionary);


// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

deSouzaMendesThompson::deSouzaMendesThompson
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    rheologyModel(U, phi),
    //viscoElasticModel(U, phi),
    //thixotropyModel(U, phi),
    G0_(this->subDict((this->type()+"Coeffs")).lookup("G0")),
    m_(this->subDict((this->type()+"Coeffs")).lookup("m")),
    tau_
    (
        IOobject
        (
            "tau",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
        //symm(nuEq_*fvc::grad(U))
    ),
    twoD_
    (
        IOobject
        (
            "twoD",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoSymm(fvc::grad(U))
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> 
deSouzaMendesThompson::divTau(const volVectorField& U) const
{
    return
    (
      fvc::div(tau_, "div(tau)")
      + fvm::laplacian(nu_, U, "fvm::laplacian(nu,U)")
      + fvc::div(nu_*dev(T(fvc::grad(U))), "div(nu*dev(T(grad(U))))")
      - fvc::div((nu_)*fvc::grad(U), "div(tau)")
    );
}


void deSouzaMendesThompson::correct()
{
    // Update rheology base model 
    rheologyModel::correct();
    structureModelPtr_->correct(phi_);

    // Retrieve structure field
    const volScalarField& lambda = structureModelPtr_->lambda();

    const dimensionedScalar lambda0 = 
        min
        (
            log(this->nu0()/this->nuInf()), 
            dimensionedScalar("lambdaMax", dimless, 1.0e6)
        );

    // Calculate transport properties
    nu_ = nuEq_*0.0 + nuInf_;
    volScalarField Gs = 
        nu_/this->nu0()*dimensionedScalar("zero", G0_.dimensions(), 0.0);
    forAll(nu_, cellI)
    {
        nu_[cellI] *= exp(lambda[cellI]);
    }
    forAll(nu_.boundaryField(), patchI)
    {
        forAll(nu_.boundaryField()[patchI], faceI)
        {
            nu_.boundaryField()[patchI][faceI] 
                *= exp(lambda.boundaryField()[patchI][faceI]);
        }
    }

    Gs += G0_;
    forAll(Gs, cellI)
    {
        Gs[cellI] *= m_.value()*exp(1.0/lambda[cellI] - 1.0/lambda0.value());
    }
    forAll(Gs.boundaryField(), patchI)
    {
        forAll(Gs.boundaryField()[patchI], faceI)
        {
            Gs.boundaryField()[patchI][faceI] 
                *= m_.value()*exp(1.0/lambda.boundaryField()[patchI][faceI] 
                          - 1.0/lambda0.value());
        }
    }

    volScalarField theta1 = (1.0 - nuInf_/nu_)*nu_/Gs;
    volScalarField theta2 = (1.0 - nuInf_/nu_)*nuInf_/Gs;
    
    // Velocity gradient tensor
    const tmp<volTensorField> tL = fvc::grad(U_);
    const volTensorField& L = tL();

    // Convected derivate terms
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    twoD_ = twoSymm(L);
    volTensorField twoDgradU = twoD_ & L;

    if(db().foundObject<solutionControl>("simpleControl"))
    {
        // Stress transport equation
        fvSymmTensorMatrix tauEqn
        (
            fvm::div(phi_, tau_)
         ==
            nuInf_*twoD_/theta2
          + nuInf_*(fvc::div(phi_, twoD_) - twoSymm(twoDgradU))
          + twoSymm(C)
          - fvm::Sp(1/theta1, tau_)
        );

        tauEqn.relax();
        tauEqn.solve();
    }
    else
    {
        // Stress transport equation
        fvSymmTensorMatrix tauEqn
        (
            fvm::ddt(tau_)
          + fvm::div(phi_, tau_)
         ==
            nuInf_*twoD_/theta2
          + nuInf_*(fvc::ddt(twoD_) + fvc::div(phi_, twoD_) - twoSymm(twoDgradU))
          + twoSymm(C)
          - fvm::Sp(1/theta1, tau_)
        );

        tauEqn.relax();
        tauEqn.solve();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

#include "rheologyModel.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rheologyModel, 0);
defineRunTimeSelectionTable(rheologyModel, dictionary)

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

//tmp<volScalarField> rheologyModel::strainRate()
//{
//    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
//}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rheologyModel::rheologyModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    //singlePhaseTransportModel(U, phi),
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    U_(U),
    phi_(phi),
    viscosityModelPtr_
    (
        viscosityModel::New("nu", *this, U, phi)
    ),
    rho_(this->lookup("rho")),
    nu0_
    (
        this->subDict("rheologyCoeffs").lookupOrDefault
        (
            "nu0",
            dimensionedScalar("nu0", dimArea/dimTime, 0.0)
        )
    ),
    nuInf_
    (
        this->subDict("rheologyCoeffs").lookupOrDefault
        (
            "nuInf",
            dimensionedScalar("nuInf", nu0_.dimensions(), 0.0)
        )
    ),
    structureModelPtr_
    (
        structureModel::New
        (
            "lambda", 
            *this, 
            U
        )
    ),
    nuEq_(viscosityModelPtr_->nu()),
    sr_
    (
        IOobject
        (
            "sr",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        viscosityModelPtr_->strainRate()
    ),
    nu_
    (
        IOobject
        (
            "nuApp",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nuEq_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> 
rheologyModel::divTau(const volVectorField& U) const
{
    return
    (
      fvm::laplacian(nu_, U, "fvm::laplacian(nu,U)")
      + fvc::div(nu_*dev(T(fvc::grad(U))), "div(nu*dev(T(grad(U))))")
    );
}

void rheologyModel::correct()
{
    sr_ = viscosityModelPtr_->strainRate();
    viscosityModelPtr_->correct();
}

bool rheologyModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

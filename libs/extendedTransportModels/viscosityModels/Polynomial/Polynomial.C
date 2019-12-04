/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Polynomial.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Polynomial, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        Polynomial,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Polynomial::calcNu() 
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    dimensionedScalar nu1(nu0_/nu0_.value());
    tmp<volScalarField> sr(strainRate());

    return
    (
        min
        (
            nu0_,
            (coeff1_ + coeff2_*tone*sr() + coeff3_*pow(tone*sr(),exp1_) 
             + coeff4_*pow(tone*sr(),exp2_))
           /(max(sr()*tone, criticalSR_*tone))*nu1
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Polynomial::Polynomial
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    PolynomialCoeffs_
    (
        viscosityProperties.subDict(typeName + "Coeffs")
    ),
    nu0_(PolynomialCoeffs_.lookup("nu0")),
    coeff1_
    (
        readScalar(PolynomialCoeffs_.lookup("coeff1"))
    ),
    coeff2_
    (
        readScalar(PolynomialCoeffs_.lookup("coeff2"))
    ),
    coeff3_
    (
        readScalar(PolynomialCoeffs_.lookup("coeff3"))
    ),
    coeff4_
    (
        readScalar(PolynomialCoeffs_.lookup("coeff4"))
    ),
    exp1_
    (
        readScalar(PolynomialCoeffs_.lookup("exp1"))
    ),
    exp2_
    (
        readScalar(PolynomialCoeffs_.lookup("exp2"))
    ),
    criticalSR_(PolynomialCoeffs_.lookup("criticalSR")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::Polynomial::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    PolynomialCoeffs_ 
        = viscosityProperties.subDict(typeName + "Coeffs");

    PolynomialCoeffs_.lookup("nu0") >> nu0_;
    PolynomialCoeffs_.lookup("coeff1") >> coeff1_;
    PolynomialCoeffs_.lookup("coeff2") >> coeff2_;
    PolynomialCoeffs_.lookup("coeff3") >> coeff3_;
    PolynomialCoeffs_.lookup("coeff4") >> coeff4_;
    PolynomialCoeffs_.lookup("exp1") >> exp1_;
    PolynomialCoeffs_.lookup("exp2") >> exp2_;
    PolynomialCoeffs_.lookup("criticalSR") >> criticalSR_;

    return true;
}


// ************************************************************************* //

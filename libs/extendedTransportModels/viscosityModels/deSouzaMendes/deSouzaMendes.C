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

#include "deSouzaMendes.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(deSouzaMendes, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        deSouzaMendes,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::deSouzaMendes::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    tmp<volScalarField> sr(strainRate());
    
    return
    (
        min
        (
            nu0_,
            (1.0 - exp(-nu0_*sr()/tauY_))
            *((tauY_ - tauYD_)*exp(-sr()/gammaYD_) 
             + tauYD_ + K_*rtone*pow(sr()*tone, n_))
            /(max(sr(),dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)))
            + nuInf_
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::deSouzaMendes::deSouzaMendes
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    deSouzaMendesCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    K_(deSouzaMendesCoeffs_.lookup("K")),
    n_(deSouzaMendesCoeffs_.lookup("n")),
    tauY_(deSouzaMendesCoeffs_.lookup("tauY")),
    tauYD_(deSouzaMendesCoeffs_.lookup("tauYD")),
    gammaYD_(deSouzaMendesCoeffs_.lookup("gammaYD")),
    nu0_(deSouzaMendesCoeffs_.lookup("nu0")),
    nuInf_(deSouzaMendesCoeffs_.lookup("nuInf")),
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

bool Foam::viscosityModels::deSouzaMendes::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    deSouzaMendesCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    deSouzaMendesCoeffs_.lookup("K") >> K_;
    deSouzaMendesCoeffs_.lookup("n") >> n_;
    deSouzaMendesCoeffs_.lookup("tauY") >> tauY_;
    deSouzaMendesCoeffs_.lookup("tauYD") >> tauYD_;
    deSouzaMendesCoeffs_.lookup("gammaYD") >> gammaYD_;
    deSouzaMendesCoeffs_.lookup("nu0") >> nu0_;
    deSouzaMendesCoeffs_.lookup("nuInf") >> nuInf_;

    return true;
}


// ************************************************************************* //

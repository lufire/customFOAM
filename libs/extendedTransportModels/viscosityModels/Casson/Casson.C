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

#include "Casson.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Casson, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        Casson,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Casson::calcNu() const
{
    dimensionedScalar vone("vone", dimArea/dimTime, 1.0);
    dimensionedScalar rvone("rvone", dimTime/dimArea, 1.0);
    tmp<volScalarField> sr(strainRate());
    
    return
    (
        min
        (
            nu0_,
            pow
            (
                pow(nuInf_*rvone, n_) 
                + pow(tau0_*(1.0-exp(-sr()/srMin_))
                      /max(sr(), srMin_)*rvone, 
                      n_),
                1.0/n_
            )*vone
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Casson::Casson
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    CassonCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nu0_(CassonCoeffs_.lookup("nu0")),
    nuInf_(CassonCoeffs_.lookup("nuInf")),
    n_(CassonCoeffs_.lookup("n")),
    tau0_(CassonCoeffs_.lookup("tau0")),
    srMin_(CassonCoeffs_.lookup("srMin")),
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

bool Foam::viscosityModels::Casson::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    CassonCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    CassonCoeffs_.lookup("nu0") >> nu0_;
    CassonCoeffs_.lookup("nuInf") >> nuInf_;
    CassonCoeffs_.lookup("n") >> n_;
    CassonCoeffs_.lookup("tau0") >> tau0_;
    CassonCoeffs_.lookup("srMin") >> srMin_;

    return true;
}


// ************************************************************************* //

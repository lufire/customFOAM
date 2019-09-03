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

#include "HerschelBulkleyKD.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(HerschelBulkleyKD, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HerschelBulkleyKD,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::HerschelBulkleyKD::calcNu() 
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    tmp<volScalarField> sr(strainRate());

    if(voidfractionHeader_.headerOk())
    {
        return
        (
            min
            (
                nu0_,
                (tau0_*(1.0-exp(-m_*tone*sr())) + k_*rtone*pow(tone*sr(), n_))
               /(max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
            )
            *pow((1.0-(1.0-voidfraction_)/maxSolidFraction_),
                -EinsteinCoeff_*maxSolidFraction_)
        );
    }
    else
    {
        return
        (
            min
            (
                nu0_,
                (tau0_*(1.0-exp(-m_*tone*sr())) + k_*rtone*pow(tone*sr(), n_))
               /(max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::HerschelBulkleyKD::HerschelBulkleyKD
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    HerschelBulkleyKDCoeffs_
    (
        viscosityProperties.subDict(typeName + "Coeffs")
    ),
    k_(HerschelBulkleyKDCoeffs_.lookup("k")),
    n_(HerschelBulkleyKDCoeffs_.lookup("n")),
    tau0_(HerschelBulkleyKDCoeffs_.lookup("tau0")),
    nu0_(HerschelBulkleyKDCoeffs_.lookup("nu0")),
    m_(HerschelBulkleyKDCoeffs_.lookup("m")),
    maxSolidFraction_
    (
        readScalar(HerschelBulkleyKDCoeffs_.lookup("maxSolidFraction"))
    ),
    EinsteinCoeff_
    (
        readScalar(HerschelBulkleyKDCoeffs_.lookup("EinsteinCoeff"))
    ),
    voidfractionName_
    (
        HerschelBulkleyKDCoeffs_.lookup("voidfractionName")
    ),
    voidfraction_
    (
        U_.db().lookupObject<volScalarField>(voidfractionName_)
    ),
    voidfractionHeader_ 
    (
        voidfractionName_,
        U_.time().timeName(),
        U_.db(),
        IOobject::NO_READ
    ),
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

bool Foam::viscosityModels::HerschelBulkleyKD::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HerschelBulkleyKDCoeffs_ 
        = viscosityProperties.subDict(typeName + "Coeffs");

    HerschelBulkleyKDCoeffs_.lookup("k") >> k_;
    HerschelBulkleyKDCoeffs_.lookup("n") >> n_;
    HerschelBulkleyKDCoeffs_.lookup("tau0") >> tau0_;
    HerschelBulkleyKDCoeffs_.lookup("nu0") >> nu0_;
    HerschelBulkleyKDCoeffs_.lookup("m") >> m_;
    HerschelBulkleyKDCoeffs_.lookup("maxSolidFraction")
        >> maxSolidFraction_;
    HerschelBulkleyKDCoeffs_.lookup("EinsteinCoeff") 
        >> EinsteinCoeff_;
    HerschelBulkleyKDCoeffs_.lookup("voidfractionName") 
        >> voidfractionName_;

    return true;
}


// ************************************************************************* //

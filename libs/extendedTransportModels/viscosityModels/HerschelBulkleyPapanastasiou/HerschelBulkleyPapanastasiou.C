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

#include "HerschelBulkleyPapanastasiou.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(HerschelBulkleyPapanastasiou, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HerschelBulkleyPapanastasiou,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::HerschelBulkleyPapanastasiou::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    tmp<volScalarField> sr(strainRate());
    
    return
    (
        min
        (
            nu0_,
            (tau0_*(1.0-exp(-m_*tone*sr())) + k_*rtone*pow(tone*sr(), n_))
           /(max(sr(), srMin_))
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::HerschelBulkleyPapanastasiou::HerschelBulkleyPapanastasiou
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    HerschelBulkleyPapanastasiouCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    k_(HerschelBulkleyPapanastasiouCoeffs_.lookup("k")),
    n_(HerschelBulkleyPapanastasiouCoeffs_.lookup("n")),
    tau0_(HerschelBulkleyPapanastasiouCoeffs_.lookup("tau0")),
    nu0_(HerschelBulkleyPapanastasiouCoeffs_.lookup("nu0")),
    m_(HerschelBulkleyPapanastasiouCoeffs_.lookup("m")),
    srMin_(HerschelBulkleyPapanastasiouCoeffs_.lookup("srMin")),
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
    /*kinVisc_
    (
	nu_
    )*/
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::HerschelBulkleyPapanastasiou::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    HerschelBulkleyPapanastasiouCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    HerschelBulkleyPapanastasiouCoeffs_.lookup("k") >> k_;
    HerschelBulkleyPapanastasiouCoeffs_.lookup("n") >> n_;
    HerschelBulkleyPapanastasiouCoeffs_.lookup("tau0") >> tau0_;
    HerschelBulkleyPapanastasiouCoeffs_.lookup("nu0") >> nu0_;
    HerschelBulkleyPapanastasiouCoeffs_.lookup("m") >> m_;
    HerschelBulkleyPapanastasiouCoeffs_.lookup("srMin") >> srMin_;

    return true;
}


// ************************************************************************* //

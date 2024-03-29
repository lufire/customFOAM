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

#include "zeroNormalGradientFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroNormalGradientFvPatchVectorField::
zeroNormalGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::zeroNormalGradientFvPatchVectorField::
zeroNormalGradientFvPatchVectorField
(
    const zeroNormalGradientFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper)
{}


Foam::zeroNormalGradientFvPatchVectorField::
zeroNormalGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=
        (
            (this->patchInternalField() & this->patch().nf())
            *this->patch().nf()
        );
}


Foam::zeroNormalGradientFvPatchVectorField::
zeroNormalGradientFvPatchVectorField
(
    const zeroNormalGradientFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf)
{}


Foam::zeroNormalGradientFvPatchVectorField::
zeroNormalGradientFvPatchVectorField
(
    const zeroNormalGradientFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroNormalGradientFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    vectorField::operator=
    (    
        (this->patch().nf() & this->patchInternalField())*this->patch().nf()
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::zeroNormalGradientFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::zeroNormalGradientFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        zeroNormalGradientFvPatchVectorField
    );
}

// ************************************************************************* //

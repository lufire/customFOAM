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

#include "electroNeutralityFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroNeutralityFvPatchScalarField::
electroNeutralityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    zCName_("zC"),
    chargeNumber_(0)
{
    const word& speciesFieldName = this->dimensionedInternalField().name();
    size_t found = speciesFieldName.find("-");
    if(found == std::string::npos)
    {
        found = speciesFieldName.find("+");
        if(found == std::string::npos)
        {
            chargeNumber_ = 0;
        }
        else
        {
            chargeNumber_ = speciesFieldName.size() - found;
        }
    }
    else
    {
        chargeNumber_ = (speciesFieldName.size() - found)*(-1);
    }
}


Foam::electroNeutralityFvPatchScalarField::
electroNeutralityFvPatchScalarField
(
    const electroNeutralityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    zCName_(ptf.zCName_),
    chargeNumber_(ptf.chargeNumber_)
{}


Foam::electroNeutralityFvPatchScalarField::
electroNeutralityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    zCName_(dict.lookupOrDefault<word>("zCField", "zC")),
    chargeNumber_(0)
{

    const word& speciesFieldName = this->dimensionedInternalField().name();
    size_t found = speciesFieldName.find("-");
    if(found == std::string::npos)
    {
        found = speciesFieldName.find("+");
        if(found == std::string::npos)
        {
            chargeNumber_ = 0;
        }
        else
        {
            chargeNumber_ = speciesFieldName.size() - found;
        }
    }
    else
    {
        chargeNumber_ = (speciesFieldName.size() - found)*(-1);
    }
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::electroNeutralityFvPatchScalarField::
electroNeutralityFvPatchScalarField
(
    const electroNeutralityFvPatchScalarField& bvppsf
)
:
    fixedValueFvPatchScalarField(bvppsf),
    zCName_(bvppsf.zCName_),
    chargeNumber_(bvppsf.chargeNumber_)
{
}

Foam::electroNeutralityFvPatchScalarField::
electroNeutralityFvPatchScalarField
(
    const electroNeutralityFvPatchScalarField& bvppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(bvppsf, iF),
    zCName_(bvppsf.zCName_),
    chargeNumber_(bvppsf.chargeNumber_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electroNeutralityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Return const access to patch values of current density field
    const volScalarField& zC =
         db().lookupObject<volScalarField>(zCName_);
    const fvPatchField<scalar>& zCp =
        patch().patchField<volScalarField, scalar>(zC);    

    // Transform vectors into the local coordinate system
              
    if (zC.dimensions() == dimMoles/dimVol)
    {
        scalarField::operator=
        (
           -zCp/chargeNumber_ 
        );
    }
    else
    {
        FatalErrorIn("electroNeutralityFvPatchScalarField::updateCoeffs()")
            << "dimensions of zC are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::electroNeutralityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "zC", "zC", zCName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
void Foam::electroNeutralityFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& psf
)
{
    fvPatchField<scalar>::operator=(patch().nf()*(patch().nf() & psf));
}
*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electroNeutralityFvPatchScalarField
    );
}

// ************************************************************************* //

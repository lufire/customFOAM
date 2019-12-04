/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "linearGradientFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearGradientFvPatchScalarField::linearGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1),
    fieldName_("default")
{}


Foam::linearGradientFvPatchScalarField::linearGradientFvPatchScalarField
(
    const linearGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1),
    fieldName_(ptf.fieldName_)

{
    patchType() = ptf.patchType();

    // Map gradient. Set unmapped values and overwrite with mapped ptf
    gradient() = 0.0;
    gradient().map(ptf.gradient(), mapper);

    // Evaluate the value field from the gradient if the internal field is valid
    if (&iF && iF.size())
    {
        scalarField::operator=
        (
            //patchInternalField() + gradient()/patch().deltaCoeffs()
            // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
            // which fails for AMI patches for some mapping operations
            patchInternalField() + gradient()*(patch().nf() & patch().delta())
        );
    }
}


Foam::linearGradientFvPatchScalarField::linearGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1),
    fieldName_(dict.lookup("fieldName"))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::linearGradientFvPatchScalarField::linearGradientFvPatchScalarField
(
    const linearGradientFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    curTimeIndex_(-1),
    fieldName_(wbppsf.fieldName_)
{}


Foam::linearGradientFvPatchScalarField::linearGradientFvPatchScalarField
(
    const linearGradientFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    curTimeIndex_(-1),
    fieldName_(wbppsf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linearGradientFvPatchScalarField::updateCoeffs
(
    const scalarField& linearCoeffp
)
{
    if (updated())
    {
        return;
    }

    curTimeIndex_ = this->db().time().timeIndex();
    const fvPatchVectorField& patchField =
       patch().lookupPatchField<volVectorField, vector>(fieldName_);

    gradient() = (linearCoeffp * patchField) & patch().nf();
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::linearGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        FatalErrorIn("linearGradientFvPatchScalarField::updateCoeffs()")
            << "updateCoeffs(const scalarField& linearCoeffp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);
    }
}


void Foam::linearGradientFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("fieldName") << fieldName_ 
        << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        linearGradientFvPatchScalarField
    );
}


// ************************************************************************* //

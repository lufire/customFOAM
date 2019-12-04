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

#include "currentSpeciesFluxFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::currentSpeciesFluxFvPatchVectorField::
currentSpeciesFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    iName_("i"),
    stoichCoeff_(1),
    electronNumber_(1)
{}


Foam::currentSpeciesFluxFvPatchVectorField::
currentSpeciesFluxFvPatchVectorField
(
    const currentSpeciesFluxFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    iName_(ptf.iName_),
    stoichCoeff_(ptf.stoichCoeff_),
    electronNumber_(ptf.electronNumber_)
{}


Foam::currentSpeciesFluxFvPatchVectorField::
currentSpeciesFluxFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    iName_(dict.lookupOrDefault<word>("currentField","i")),
    stoichCoeff_(readLabel(dict.lookup("stoichCoeff"))),
    electronNumber_(readLabel(dict.lookup("electronNumber")))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::currentSpeciesFluxFvPatchVectorField::
currentSpeciesFluxFvPatchVectorField
(
    const currentSpeciesFluxFvPatchVectorField& bvcpvf
)
:
    fixedValueFvPatchVectorField(bvcpvf),
    iName_(bvcpvf.iName_),
    stoichCoeff_(bvcpvf.stoichCoeff_),
    electronNumber_(bvcpvf.electronNumber_)
{}


Foam::currentSpeciesFluxFvPatchVectorField::
currentSpeciesFluxFvPatchVectorField
(
    const currentSpeciesFluxFvPatchVectorField& bvcpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(bvcpvf, iF),
    iName_(bvcpvf.iName_),
    stoichCoeff_(bvcpvf.stoichCoeff_),
    electronNumber_(bvcpvf.electronNumber_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::currentSpeciesFluxFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Physical constants
    //const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();

    //const dictionary& thermoDict = db().lookupObject<IOdictionary>
    //(
    //    "thermophysicalProperties"
    //);   
    // Return const access to internal cell values adjacent to
    // patch of electric field
    const volVectorField& i =
         db().lookupObject<volVectorField>(iName_); 

    const fvPatchField<vector>& ip =
        patch().patchField<volVectorField, vector>(i);
 
    //const tmp<vectorField> tipI =
    //      patch().patchInternalField(i.internalField());
    //      
    //const vectorField& ipI = tipI();
    
    //const tmp<vectorField>& tn = patch().nf();
    //const vectorField& n = tn();
//    const Field<scalar>& magS = patch().magSf();
                  
    if (i.dimensions() == dimCurrent/dimArea)
    {
        vectorField::operator=
        (    
            -ip*stoichCoeff_/(electronNumber_*F)
        );
    }
    else
    {
        FatalErrorIn("currentSpeciesFluxFvPatchVectorField::updateCoeffs()")
            << "dimensions of i are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::currentSpeciesFluxFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>
    (
        os, 
        "currentField", 
        "i", 
        iName_
    );
    writeEntry(os, "stoichCoeff", stoichCoeff_);
    writeEntry(os, "electronNumber", electronNumber_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::currentSpeciesFluxFvPatchVectorField::operator=
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
        currentSpeciesFluxFvPatchVectorField
    );
}

// ************************************************************************* //

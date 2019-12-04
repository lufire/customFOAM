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

#include "shearStressSlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shearStressSlipVelocityFvPatchVectorField::
shearStressSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    nuName_("nu"),
    factor_(1.0),
    exponent_(1.0),
    rho_(0.0)
{
    const IOdictionary transDict 
    (
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    rho_ = transDict.lookup("rho");
}


Foam::shearStressSlipVelocityFvPatchVectorField::
shearStressSlipVelocityFvPatchVectorField
(
    const shearStressSlipVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    nuName_(pvf.nuName_),
    factor_(pvf.factor_),
    exponent_(pvf.exponent_),
    rho_(pvf.rho_)
{}


Foam::shearStressSlipVelocityFvPatchVectorField::
shearStressSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    nuName_(dict.lookupOrDefault<word>("nuName","nu")),
    factor_(readScalar(dict.lookup("factor"))),
    exponent_(readScalar(dict.lookup("exponent"))),
    rho_
    (
        dict.lookupOrDefault<dimensionedScalar>
        (
            "rho",
            dimensionedScalar("zero", dimMass/dimVol, 0.0)
        )
    )
{
    const IOdictionary transDict 
    (
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    rho_ = transDict.lookup("rho");
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::shearStressSlipVelocityFvPatchVectorField::
shearStressSlipVelocityFvPatchVectorField
(
    const shearStressSlipVelocityFvPatchVectorField& pvf 
)
:
    fixedValueFvPatchVectorField(pvf),
    nuName_(pvf.nuName_),
    factor_(pvf.factor_),
    exponent_(pvf.exponent_),
    rho_(pvf.rho_)
{}


Foam::shearStressSlipVelocityFvPatchVectorField::
shearStressSlipVelocityFvPatchVectorField
(
    const shearStressSlipVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    nuName_(pvf.nuName_),
    factor_(pvf.factor_),
    exponent_(pvf.exponent_),
    rho_(pvf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shearStressSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const volVectorField& U =
         db().lookupObject<volVectorField>("U");

    const volVectorField& Uprev = U.prevIter();

    const fvPatchField<vector>& Uprevp =
        patch().patchField<volVectorField, vector>(Uprev);
    
    const tmp<vectorField> UprevpI =
          patch().patchInternalField(Uprev.internalField());
    
    const volScalarField& nu = 
        db().lookupObject<volScalarField>(nuName_); 
    const volScalarField& nuPrev = nu.prevIter();
    const scalarField& nuPrevp =
        patch().patchField<volScalarField, scalar>(nuPrev);
    
    tmp<vectorField> tn = patch().nf();
    const vectorField& n = tn();

    const vectorField tUprevp = Uprevp - n*(n & Uprevp);
    const vectorField tUprevpI = UprevpI() - n*(n & UprevpI());
    const vectorField nUgradPrev = (tUprevp - tUprevpI)*patch().deltaCoeffs();
    
    const scalarField d = factor_*pow(max(mag(nUgradPrev),VSMALL), (exponent_-1.0))
        *pow(max(nuPrevp*rho_.value(),VSMALL), exponent_)*patch().deltaCoeffs();
    
    vectorField Up = d/(d+1.0)*tUprevpI;

    operator==(Up);


    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::shearStressSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>
    (
        os, 
        "nuName", 
        "nu", 
        nuName_
    );
    writeEntry(os, "factor", factor_);
    writeEntry(os, "exponent", exponent_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//void Foam::shearStressSlipVelocityFvPatchVectorField::operator=
//(
//    const fvPatchField<vector>& pvf
//)
//{
//    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        shearStressSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

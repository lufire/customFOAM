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

#include "diffusionCurrentFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionCurrentFvPatchVectorField::
diffusionCurrentFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    electrodeName_("none"),
    phiEName_("phiE"),
    TName_("T"),
    oxidantName_("oxidant"),
    reductantName_("reductant"),
    CRefOxidant_(1.0),
    CRefReductant_(1.0),
    gammaOxidant_(1.0),
    gammaReductant_(1.0),
    jEx_(0.0),
    equilibriumPotential_(0.0),
    alphaA_(1.0),
    alphaC_(1.0)
{}


Foam::diffusionCurrentFvPatchVectorField::
diffusionCurrentFvPatchVectorField
(
    const diffusionCurrentFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    electrodeName_(ptf.electrodeName_),
    phiEName_(ptf.phiEName_),
    TName_(ptf.TName_),
    oxidantName_(ptf.oxidantName_),
    reductantName_(ptf.reductantName_),
    CRefOxidant_(ptf.CRefOxidant_),
    CRefReductant_(ptf.CRefReductant_),
    gammaOxidant_(ptf.gammaOxidant_),
    gammaReductant_(ptf.gammaReductant_),
    jEx_(ptf.jEx_),
    equilibriumPotential_(ptf.equilibriumPotential_),
    electrodePotential_(ptf.electrodePotential_),
    alphaA_(ptf.alphaA_),
    alphaC_(ptf.alphaC_)
{}


Foam::diffusionCurrentFvPatchVectorField::
diffusionCurrentFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    electrodeName_(dict.lookupOrDefault<word>("electrodeName", "none")),
    phiEName_(dict.lookupOrDefault<word>("electricPotentialField", "phiE")),
    TName_(dict.lookupOrDefault<word>("temperatureField","T")),
    oxidantName_(dict.lookup("oxidant")),
    reductantName_(dict.lookup("reductant")),
    CRefOxidant_(dict.lookupOrDefault<scalar>("CRefOxidant", 1.0)),
    CRefReductant_(dict.lookupOrDefault<scalar>("CRefReductant", 1.0)),
    gammaOxidant_(dict.lookupOrDefault<scalar>("gammaOxidant", 1.0)),
    gammaReductant_(dict.lookupOrDefault<scalar>("gammaReductant", 1.0)),
    jEx_(readScalar(dict.lookup("exchangeCurrentDensity"))),
    equilibriumPotential_(dict.lookupOrDefault<scalar>("equilibriumPotential", 0.0)),
    electrodePotential_(dict.lookupOrDefault<scalar>("electrodePotential", 0.0)),
    alphaA_(dict.lookupOrDefault<scalar>("alphaA", 1.0)),
    alphaC_(dict.lookupOrDefault<scalar>("alphaC", 1.0))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::diffusionCurrentFvPatchVectorField::
diffusionCurrentFvPatchVectorField
(
    const diffusionCurrentFvPatchVectorField& bvcpvf
)
:
    fixedValueFvPatchVectorField(bvcpvf),
    electrodeName_(bvcpvf.electrodeName_),
    phiEName_(bvcpvf.phiEName_),
    TName_(bvcpvf.TName_),
    oxidantName_(bvcpvf.oxidantName_),
    reductantName_(bvcpvf.reductantName_),
    CRefOxidant_(bvcpvf.CRefOxidant_),
    CRefReductant_(bvcpvf.CRefReductant_),
    gammaOxidant_(bvcpvf.gammaOxidant_),
    gammaReductant_(bvcpvf.gammaReductant_),
    jEx_(bvcpvf.jEx_),
    equilibriumPotential_(bvcpvf.equilibriumPotential_),
    electrodePotential_(bvcpvf.electrodePotential_),
    alphaA_(bvcpvf.alphaA_),
    alphaC_(bvcpvf.alphaC_)
{}


Foam::diffusionCurrentFvPatchVectorField::
diffusionCurrentFvPatchVectorField
(
    const diffusionCurrentFvPatchVectorField& bvcpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(bvcpvf, iF),
    electrodeName_(bvcpvf.electrodeName_),
    phiEName_(bvcpvf.phiEName_),
    TName_(bvcpvf.TName_),
    oxidantName_(bvcpvf.oxidantName_),
    reductantName_(bvcpvf.reductantName_),
    CRefOxidant_(bvcpvf.CRefOxidant_),
    CRefReductant_(bvcpvf.CRefReductant_),
    gammaOxidant_(bvcpvf.gammaOxidant_),
    gammaReductant_(bvcpvf.gammaReductant_),
    jEx_(bvcpvf.jEx_),
    equilibriumPotential_(bvcpvf.equilibriumPotential_),
    electrodePotential_(bvcpvf.electrodePotential_),
    alphaA_(bvcpvf.alphaA_),
    alphaC_(bvcpvf.alphaC_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diffusionCurrentFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();

//     // Return const access to patch values of electric field
//     const surfaceScalarField& phiEs =
//         db().lookupObject<surfaceScalarField>(phiEName_);
// 
//     const fvsPatchField<scalar>& phiEp =
//         patch().patchField<surfaceScalarField, scalar>(phiEs);

    const scalarField& deltaCoeffs = patch().deltaCoeffs();


    // Return const access to internal cell values adjacent to
    // patch of electric field
    const volScalarField& phiE =
         db().lookupObject<volScalarField>(phiEName_); 

    const fvPatchField<scalar>& phiEp =
        patch().patchField<volScalarField, scalar>(phiE);
 
    const tmp<scalarField> tphiEpI =
          patch().patchInternalField(phiE.internalField());
          
    const scalarField& phiEpI = tphiEpI();
    
     const volScalarField& sigma =
        db().lookupObject<volScalarField>("sigma");
   
    const tmp<vectorField>& tn = patch().nf();
    const vectorField& n = tn();
//    const Field<scalar>& magS = patch().magSf();
                  
    if (phiE.dimensions() == dimMass*dimArea/pow3(dimTime)/dimCurrent)
    {
        vectorField::operator=
        (    
            (-sigma*(phiEpI-phiEp)*deltaCoeffs)*n*(-1.0)
        );
    }
    else
    {
        FatalErrorIn("diffusionCurrentFvPatchVectorField::updateCoeffs()")
            << "dimensions of phiE are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::diffusionCurrentFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phiE", "phiE", phiEName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::diffusionCurrentFvPatchVectorField::operator=
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
        diffusionCurrentFvPatchVectorField
    );
}

// ************************************************************************* //

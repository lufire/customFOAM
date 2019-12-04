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

#include "fixedCurrentConcentrationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"
//#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedCurrentConcentrationFvPatchScalarField::
fixedCurrentConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    electrodeName_("none"),
    jName_("j"),
    TName_("T"),
    oxidantName_("oxidant"),
    reductantName_("reductant"),
    CRefOxidant_(1.0),
    CRefReductant_(1.0),
    gammaOxidant_(1.0),
    gammaReductant_(1.0),
    jEx_(0.0),
    equilibriumPotential_(0.0),
    valenceNumber_(1.0),
    alphaC_(1.0)
{}


Foam::fixedCurrentConcentrationFvPatchScalarField::
fixedCurrentConcentrationFvPatchScalarField
(
    const fixedCurrentConcentrationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    electrodeName_(ptf.electrodeName_),
    jName_(ptf.jName_),
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
    valenceNumber_(ptf.valenceNumber_),
    alphaC_(ptf.alphaC_)
{}


Foam::fixedCurrentConcentrationFvPatchScalarField::
fixedCurrentConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    electrodeName_(dict.lookupOrDefault<word>("electrodeName", "none")),
    jName_(dict.lookupOrDefault<word>("currentDensityField", "j")),
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
    valenceNumber_(dict.lookupOrDefault<scalar>("valenceNumber", 1.0)),
    alphaC_(dict.lookupOrDefault<scalar>("alphaC", 1.0))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::fixedCurrentConcentrationFvPatchScalarField::
fixedCurrentConcentrationFvPatchScalarField
(
    const fixedCurrentConcentrationFvPatchScalarField& bvppsf
)
:
    fixedValueFvPatchScalarField(bvppsf),
    electrodeName_(bvppsf.electrodeName_),
    jName_(bvppsf.jName_),
    TName_(bvppsf.TName_),
    oxidantName_(bvppsf.oxidantName_),
    reductantName_(bvppsf.reductantName_),
    CRefOxidant_(bvppsf.CRefOxidant_),
    CRefReductant_(bvppsf.CRefReductant_),
    gammaOxidant_(bvppsf.gammaOxidant_),
    gammaReductant_(bvppsf.gammaReductant_),
    jEx_(bvppsf.jEx_),
    equilibriumPotential_(bvppsf.equilibriumPotential_),
    electrodePotential_(bvppsf.electrodePotential_),
    valenceNumber_(bvppsf.valenceNumber_),
    alphaC_(bvppsf.alphaC_)
{}


Foam::fixedCurrentConcentrationFvPatchScalarField::
fixedCurrentConcentrationFvPatchScalarField
(
    const fixedCurrentConcentrationFvPatchScalarField& bvppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(bvppsf, iF),
    electrodeName_(bvppsf.electrodeName_),
    jName_(bvppsf.jName_),
    TName_(bvppsf.TName_),
    oxidantName_(bvppsf.oxidantName_),
    reductantName_(bvppsf.reductantName_),
    CRefOxidant_(bvppsf.CRefOxidant_),
    CRefReductant_(bvppsf.CRefReductant_),
    gammaOxidant_(bvppsf.gammaOxidant_),
    gammaReductant_(bvppsf.gammaReductant_),
    jEx_(bvppsf.jEx_),
    equilibriumPotential_(bvppsf.equilibriumPotential_),
    electrodePotential_(bvppsf.electrodePotential_),
    valenceNumber_(bvppsf.valenceNumber_),
    alphaC_(bvppsf.alphaC_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::fixedCurrentConcentrationFvPatchScalarField::bvNewtonMethodFunction
(
    const scalar& overPotential,
    const label& faceI
)
{

}

void Foam::fixedCurrentConcentrationFvPatchScalarField::bvNewtonMethodDerivative
(
    const scalar& overPotential
)
{

}
*/
    

void Foam::fixedCurrentConcentrationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
   
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();
    
    const scalarField& deltaCoeffs = this->patch().deltaCoeffs();
    
    const dictionary& thermoDict = db().lookupObject<IOdictionary>
    (
        "thermophysicalProperties"
    );   

    // Return const access to internal cell values adjacent to
    // patch of electric field
    
    const volScalarField& C =
         db().lookupObject<volScalarField>
         (this->dimensionedInternalField().name()); 

    const fvPatchField<scalar>& Cp =
        patch().patchField<volScalarField, scalar>(C);
 
    const tmp<scalarField> tCpI =
          patch().patchInternalField(C.internalField());
          

    const scalarField& CpI = tCpI();
    

    // Return const access to patch values of current density field
    const volScalarField& phiE =
         db().lookupObject<volScalarField>("phiE");

    const fvPatchField<scalar>& phiEp =
        patch().patchField<volScalarField, vector>(phiE);

    const tmp<scalarField> phiEpI =
          patch().patchInternalField(phiE.internalField());

    // Return const access to patch values of current density field
    const volVectorField& j =
         db().lookupObject<volVectorField>(jName_);

    const fvPatchField<vector>& jp =
        patch().patchField<volVectorField, vector>(j);
    

    const tmp<vectorField>& tn = patch().nf();
    vectorField n = tn();
//    const Field<scalar>& magS = patch().magSf();
    
    if (j.dimensions() == dimCurrent/dimArea)
    {
        if(this->dimensionedInternalField().name() == "C_"+oxidantName_)
        {
            const volScalarField& D =
                db().lookupObject<volScalarField>("D_"+oxidantName_);

            const fvPatchField<scalar>& Dp =
                patch().patchField<volScalarField, scalar>(D);

            //Cp = CpI + (jp&n)/(Dp*F*valenceNumber_*this->patch().deltaCoeffs())*(-1.0);
            scalarField::operator=
            (
                this->patchInternalField() + (jp&n)
                   /(Dp*F*valenceNumber_*deltaCoeffs)*(-1.0)
            );            
            Info << "Name of Species at boundary patch: " << this->dimensionedInternalField().name() << endl;
            Info << "Patch-cell distance: " << 1.0/deltaCoeffs[0] << endl;
            Info << "Diffusion coefficient: " << Dp[0] << endl;
            Info << "Patch concentration: " <<  CpI[0] + (jp[0]&n[0])/(Dp[0]*F*valenceNumber_*deltaCoeffs[0])*(-1.0) << endl;
            Info << "Patch internal concentration: " << CpI[0] << endl;
            Info << "Patch current density: " << jp[0] << endl;
                      
//            operator== this->patchInternalField()  
//                + (jp&n)
//                /(D*F*valenceNumber_*this->patch().deltaCoeffs())*(-1.0);
        }
        else if(this->dimensionedInternalField().name() == "C_"+reductantName_)
        {
            const volScalarField& D =
                db().lookupObject<volScalarField>("D_"+reductantName_);
            
            const fvPatchField<scalar>& Dp =
                patch().patchField<volScalarField, scalar>(D);

            //Cp = CpI + (jp&n)/(Dp*F*valenceNumber_*this->patch().deltaCoeffs())*(-1.0);

            scalarField::operator=
            (
                this->patchInternalField() - (jp&n)
                    /(Dp*F*valenceNumber_*deltaCoeffs)*(-1.0)
            );
            Info << "Name of Species at boundary patch: " << this->dimensionedInternalField().name() << endl;
            Info << "Patch-cell distance: " << 1.0/deltaCoeffs[0] << endl;
            Info << "Diffusion coefficient: " << Dp[0] << endl;
            Info << "Patch concentration: " <<  CpI[0] + (jp[0]&n[0])/(Dp[0]*F*valenceNumber_*deltaCoeffs[0])*(-1.0) << endl;
            Info << "Patch internal concentration: " << CpI[0] << endl;
            Info << "Patch current density: " << jp[0] << endl;
//            operator== this->patchInternalField()  
//                - (jp&n)
//                /(D*F*valenceNumber_*this->patch().deltaCoeffs())*(-1.0);
        }
        else if(this->dimensionedInternalField().name() == "C_K+")
        {
            
            const volScalarField& D =
                db().lookupObject<volScalarField>("D_"+reductantName_);
            
            const fvPatchField<scalar>& Dp =
                patch().patchField<volScalarField, scalar>(D);

            scalarField::operator=
            (
                this->patchInternalField()
            );
            Info << "Name of Species at boundary patch: " << this->dimensionedInternalField().name() << endl;
            Info << "Patch-cell distance: " << 1.0/deltaCoeffs[0] << endl;
            Info << "Diffusion coefficient: " << Dp[0] << endl;
            Info << "Patch concentration: " <<  this->patchInternalField() << endl;
            Info << "Patch internal concentration: " << CpI[0] << endl;
            Info << "Patch current density: " << jp[0] << endl;
//            FatalErrorIn("fixedCurrentConcentrationFvPatchScalarField::updateCoeffs()")
//                << "Name of concentration field does"
//                << "not correspond to specified oxidant or reductant name"
//                << "\n    on patch " << this->patch().name()
//                << " of field " << this->dimensionedInternalField().name()
//                << " in file " << this->dimensionedInternalField().objectPath()
//                << exit(FatalError);
        }
        else
        {
            FatalErrorIn("fixedCurrentConcentrationFvPatchScalarField::updateCoeffs()")
                << "Name of concentration field does"
                << "not correspond to specified oxidant or reductant name"
                << "\n    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("fixedCurrentConcentrationFvPatchScalarField::updateCoeffs()")
            << "dimensions of j are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedCurrentConcentrationFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "j", "j", jName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
void Foam::fixedCurrentConcentrationFvPatchScalarField::operator=
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
        fixedCurrentConcentrationFvPatchScalarField
    );
}

// ************************************************************************* //

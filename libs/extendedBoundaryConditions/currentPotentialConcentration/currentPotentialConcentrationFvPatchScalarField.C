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

#include "currentPotentialConcentrationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"
//#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Foam::currentPotentialConcentrationFvPatchScalarField::
//currentPotentialConcentrationFvPatchScalarField
//(
//    const fvPatch& p,
//    const DimensionedField<scalar, volMesh>& iF
//)
//:
//    fixedValueFvPatchScalarField(p, iF),
//    currentDensityName_("i"),
//    phiEName_("phiE"),
//    TName_("T"),
//    chargeNumber_(0),
//    stoichCoeff_(0),
//    nElectrons_(1),
//    rotationMatrix_(patch().nf()->size())
//{
//    const tmp<vectorField>& tn = patch().nf();
//    const vectorField& n = tn();
//    const scalarField theta = acos(vector(1.0,0,0)&n);
//    const vectorField axis = vector(1.0,0,0)^n;
//    const scalarField a = cos(theta*0.5);
//    const scalarField b = (axis&vector(1.0,0,0))*sin(theta*0.5);
//    const scalarField c = (axis&vector(0,1.0,0))*sin(theta*0.5);
//    const scalarField d = (axis&vector(0,0,1.0))*sin(theta*0.5);
//    //rotationMatrix_(n.size()); 
//    forAll(rotationMatrix_, i)
//    {
//        if (theta[i] == 0.0)
//        {    
//            rotationMatrix_[i] = tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0); 
//        }
//        else if (theta[i] == constant::mathematical::pi)
//        {
//            rotationMatrix_[i] = tensor(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0); 
//        }
//        else
//        {
//            rotationMatrix_[i] = tensor
//            (
//                a[i]*a[i]+b[i]*b[i]-c[i]*c[i]-d[i]*d[i], 2.0*(b[i]*c[i]-a[i]*d[i]), 2.0*(b[i]*d[i]+a[i]*c[i]),
//                2.0*(b[i]*c[i]+a[i]*d[i]), a[i]*a[i]+c[i]*c[i]-b[i]*b[i]-d[i]*d[i], 2.0*(c[i]*d[i]-a[i]*b[i]),
//                2.0*(b[i]*d[i]-a[i]*c[i]), 2.0*(c[i]*d[i]+a[i]*b[i]), a[i]*a[i]+d[i]*d[i]-b[i]*b[i]-c[i]*d[i]
//            );
//        }
//    }
//}


Foam::currentPotentialConcentrationFvPatchScalarField::
currentPotentialConcentrationFvPatchScalarField
(
    const currentPotentialConcentrationFvPatchScalarField& cpcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(cpcpsf, p, iF, mapper),
    currentDensityName_(cpcpsf.currentDensityName_),
    phiEName_(cpcpsf.phiEName_),
    TName_(cpcpsf.TName_),
    chargeNumber_(cpcpsf.chargeNumber_),
    stoichCoeff_(cpcpsf.stoichCoeff_().clone().ptr()),
    nElectrons_(cpcpsf.nElectrons_().clone().ptr()),
    rotationMatrix_(cpcpsf.rotationMatrix_)
{}


//Foam::currentPotentialConcentrationFvPatchScalarField::
//currentPotentialConcentrationFvPatchScalarField
//(
//    const fvPatch& p,
//    const DimensionedField<scalar, volMesh>& iF,
//    const dictionary& dict
//)
//:
//    fixedValueFvPatchScalarField(p, iF),
//    currentDensityName_(dict.lookupOrDefault<word>("currentDensityField", "i")),
//    phiEName_(dict.lookupOrDefault<word>("electricPotentialField", "phiE")),
//    TName_(dict.lookupOrDefault<word>("temperatureField","T")),
//    chargeNumber_(dict.lookupOrDefault<label>("chargeNumber", 0)),
//    stoichCoeff_(dict.lookupOrDefault<label>("stoichCoeff", 0)),
//    nElectrons_(dict.lookupOrDefault<label>("nElectrons", 1)),
//    rotationMatrix_(patch().nf()->size())
//{
//    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
//    const tmp<vectorField>& tn = patch().nf();
//    const vectorField& n = tn();
//    const scalarField theta = acos(vector(1.0,0,0)&n);
//    const vectorField axis = vector(1.0,0,0)^n;
//    const scalarField a = cos(theta*0.5);
//    const scalarField b = (axis&vector(1.0,0,0))*sin(theta*0.5);
//    const scalarField c = (axis&vector(0,1.0,0))*sin(theta*0.5);
//    const scalarField d = (axis&vector(0,0,1.0))*sin(theta*0.5);
//    //rotationMatrix_(n.size()); 
//    forAll(rotationMatrix_, i)
//    {
//        if (theta[i] == 0.0)
//        {    
//            rotationMatrix_[i] = tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0); 
//        }
//        else if (theta[i] == constant::mathematical::pi)
//        {
//            rotationMatrix_[i] = tensor(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0); 
//        }
//        else
//        {
//            rotationMatrix_[i] = tensor
//            (
//                a[i]*a[i]+b[i]*b[i]-c[i]*c[i]-d[i]*d[i], 2.0*(b[i]*c[i]-a[i]*d[i]), 2.0*(b[i]*d[i]+a[i]*c[i]),
//                2.0*(b[i]*c[i]+a[i]*d[i]), a[i]*a[i]+c[i]*c[i]-b[i]*b[i]-d[i]*d[i], 2.0*(c[i]*d[i]-a[i]*b[i]),
//                2.0*(b[i]*d[i]-a[i]*c[i]), 2.0*(c[i]*d[i]+a[i]*b[i]), a[i]*a[i]+d[i]*d[i]-b[i]*b[i]-c[i]*d[i]
//            );
//        }
//    }
//}

Foam::currentPotentialConcentrationFvPatchScalarField::
currentPotentialConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    currentDensityName_(dict.lookup("currentDensityField")),
    phiEName_(dict.lookup("electricPotentialField")),
    TName_(dict.lookup("temperatureField")),
    chargeNumber_(0),
    stoichCoeff_(),
    nElectrons_(),
    rotationMatrix_(patch().nf()->size())
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

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

    stoichCoeff_ = DataEntry<label>::New("stoichCoeff", dict);
    nElectrons_ = DataEntry<label>::New("electronNumber", dict);

    const tmp<vectorField>& tn = patch().nf();
    const vectorField& n = tn();
    const scalarField theta = acos(vector(1.0,0,0)&n);
    const vectorField axis = vector(1.0,0,0)^n;
    const scalarField a = cos(theta*0.5);
    const scalarField b = (axis&vector(1.0,0,0))*sin(theta*0.5);
    const scalarField c = (axis&vector(0,1.0,0))*sin(theta*0.5);
    const scalarField d = (axis&vector(0,0,1.0))*sin(theta*0.5);
    //rotationMatrix_(n.size()); 
    forAll(rotationMatrix_, i)
    {
        if (theta[i] == 0.0)
        {    
            rotationMatrix_[i] = tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0); 
        }
        else if (theta[i] == constant::mathematical::pi)
        {
            rotationMatrix_[i] = tensor(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0); 
        }
        else
        {
            rotationMatrix_[i] = tensor
            (
                a[i]*a[i]+b[i]*b[i]-c[i]*c[i]-d[i]*d[i], 2.0*(b[i]*c[i]-a[i]*d[i]), 2.0*(b[i]*d[i]+a[i]*c[i]),
                2.0*(b[i]*c[i]+a[i]*d[i]), a[i]*a[i]+c[i]*c[i]-b[i]*b[i]-d[i]*d[i], 2.0*(c[i]*d[i]-a[i]*b[i]),
                2.0*(b[i]*d[i]-a[i]*c[i]), 2.0*(c[i]*d[i]+a[i]*b[i]), a[i]*a[i]+d[i]*d[i]-b[i]*b[i]-c[i]*d[i]
            );
        }
    }
}

Foam::currentPotentialConcentrationFvPatchScalarField::
currentPotentialConcentrationFvPatchScalarField
(
    const currentPotentialConcentrationFvPatchScalarField& cpcpsf
)
:
    fixedValueFvPatchScalarField(cpcpsf),
    currentDensityName_(cpcpsf.currentDensityName_),
    phiEName_(cpcpsf.phiEName_),
    TName_(cpcpsf.TName_),
    chargeNumber_(cpcpsf.chargeNumber_),
    stoichCoeff_(cpcpsf.stoichCoeff_().clone().ptr()),
    nElectrons_(cpcpsf.nElectrons_().clone().ptr()),
    rotationMatrix_(cpcpsf.rotationMatrix_)
{}


Foam::currentPotentialConcentrationFvPatchScalarField::
currentPotentialConcentrationFvPatchScalarField
(
    const currentPotentialConcentrationFvPatchScalarField& cpcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(cpcpsf, iF),
    currentDensityName_(cpcpsf.currentDensityName_),
    phiEName_(cpcpsf.phiEName_),
    TName_(cpcpsf.TName_),
    chargeNumber_(cpcpsf.chargeNumber_),
    stoichCoeff_(cpcpsf.stoichCoeff_().clone().ptr()),
    nElectrons_(cpcpsf.nElectrons_().clone().ptr()),
    rotationMatrix_(cpcpsf.rotationMatrix_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::currentPotentialConcentrationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();
    
    const scalarField delta = 1.0/this->patch().deltaCoeffs();

    const scalar t = db().time().timeOutputValue();
    const scalar stoichCoeff = stoichCoeff_->value(t);
    const scalar nElectrons = nElectrons_->value(t);
    
    word speciesName = this->dimensionedInternalField().name(); 
    speciesName = speciesName.erase(0,2);
 
    //const volScalarField& C =
    //     db().lookupObject<volScalarField>
    //     (this->dimensionedInternalField().name()); 

    //const tmp<scalarField> tCpI =
    //      patch().patchInternalField(C.internalField());

    //const scalarField& CpI = tCpI();
    

    // Return const access to patch values of the electric potential field
    const volScalarField& phiE =
         db().lookupObject<volScalarField>(phiEName_);

    const fvPatchField<scalar>& phiEp =
        patch().patchField<volScalarField, scalar>(phiE);

    const tmp<scalarField> tphiEpI =
          patch().patchInternalField(phiE.internalField());
    const scalarField& phiEpI = tphiEpI();

    // Return const access to patch values of the temperature field
    const volScalarField& T =
         db().lookupObject<volScalarField>(TName_);

    const fvPatchField<scalar>& Tp =
        patch().patchField<volScalarField, scalar>(T);

    // Return const access to patch values of the current density field
    const volVectorField& currentDensity =
         db().lookupObject<volVectorField>(currentDensityName_);

    const fvPatchField<vector>& currentDensityPatch =
        patch().patchField<volVectorField, vector>(currentDensity);
    
    //const tmp<vectorField> tjpI =
    //      patch().patchInternalField(j.internalField());
    //const vectorField& jpI = tjpI();
    //const vectorField jMean = (jp+jpI)*0.5;

    //const tmp<vectorField>& tn = patch().nf();
    //const vectorField& n = tn();

    vectorField currentDensityLocal = currentDensityPatch; 

    forAll(rotationMatrix_, i)
    {
        currentDensityLocal[i] = rotationMatrix_[i]&currentDensityPatch[i];
    }
     
    const scalarField currentDensityNormal = currentDensityLocal&vector(1.0,0,0);
    
    if (currentDensity.dimensions() == dimCurrent/dimArea)
    {
        const volScalarField& D =
            db().lookupObject<volScalarField>("D_" + speciesName);

        const fvPatchField<scalar>& Dp =
            patch().patchField<volScalarField, scalar>(D);

        const tmp<scalarField> tmup = Dp/(R*Tp);
        const scalarField& mup = tmup(); 

        const tmp<scalarField> tCpI = this->patchInternalField(); 
        const scalarField& CpI = tCpI();
        
        const tmp<scalarField> tCp = (this->patchInternalField()*Dp/delta + double(stoichCoeff)
                *currentDensityNormal/(double(nElectrons)*F))
                /(Dp/delta + double(chargeNumber_)*mup*F*(phiEp-phiEpI)/delta);
        const scalarField& Cp = tCp();

        scalarField::operator=
        (
            (this->patchInternalField()*Dp/delta - double(stoichCoeff)
                *currentDensityNormal/(double(nElectrons)*F))
                /(Dp/delta + double(chargeNumber_)*mup*F*(phiEp-phiEpI)/delta)
        );            
        Info << "Patch name: " << this->patch().name() << endl;
        Info << "Name of Species at boundary patch: " << speciesName << endl;
        Info << "Patch-cell distance: " << delta[0] << endl;
        Info << "Mobility: " << mup[0] << endl;
        Info << "Stoich coefficient: " << stoichCoeff << endl;
        Info << "Number of electrons: " << nElectrons << endl;
        Info << "Charge number: " << chargeNumber_ << endl;
        Info << "Diffusion coefficient: " << Dp[0] << endl;
        Info << "Patch internal concentration: " << CpI[0] << endl;
        Info << "Patch concentration: " << Cp[0] << endl;
        Info << "Patch electric potential: " << phiEp[0] << endl;
        Info << "Patch internal electric potential: " << phiEpI[0] << endl;
        Info << "Patch current density: " << currentDensityNormal[0] << endl;
    }
    else
    {
        FatalErrorIn("currentPotentialConcentrationFvPatchScalarField::updateCoeffs()")
            << "dimensions of current density field are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::currentPotentialConcentrationFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>
    (
        os, 
        "currentDensityField", 
        "currentDensityField", 
        currentDensityName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "electricPotentialField", 
        "electricPotentialField", 
        phiEName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "temperatureField", 
        "temperatureField", 
        TName_
    );
    stoichCoeff_->writeData(os);
    nElectrons_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//void Foam::currentPotentialConcentrationFvPatchScalarField::operator=
//(
//    const fvPatchField<scalar>& psf
//)
//{
//    fvPatchField<scalar>::operator=(patch().nf()*(patch().nf() & psf));
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        currentPotentialConcentrationFvPatchScalarField
    );
}

// ************************************************************************* //

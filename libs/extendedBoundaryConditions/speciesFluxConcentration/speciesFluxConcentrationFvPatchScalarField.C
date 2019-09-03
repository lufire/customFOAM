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

#include "speciesFluxConcentrationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"
//#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesFluxConcentrationFvPatchScalarField::
speciesFluxConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fluxFieldName_("N_"),
    phiEName_("phiE"),
    TName_("T"),
    UName_("U"),
    chargeNumber_(0),
    rotationMatrix_(patch().nf()->size())
{
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


Foam::speciesFluxConcentrationFvPatchScalarField::
speciesFluxConcentrationFvPatchScalarField
(
    const speciesFluxConcentrationFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    fluxFieldName_(psf.fluxFieldName_),
    phiEName_(psf.phiEName_),
    TName_(psf.TName_),
    UName_(psf.UName_),
    chargeNumber_(psf.chargeNumber_),
    rotationMatrix_(psf.rotationMatrix_)
{}

Foam::speciesFluxConcentrationFvPatchScalarField::
speciesFluxConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    fluxFieldName_(dict.lookupOrDefault<word>("fluxField","N_")),
    phiEName_(dict.lookupOrDefault<word>("electricPotentialField","phiE")),
    TName_(dict.lookupOrDefault<word>("temperatureField","T")),
    UName_(dict.lookupOrDefault<word>("velocityField","U")),
    chargeNumber_(0),
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

    const tmp<vectorField>& tn = patch().nf();
    const vectorField& n = tn();
    const scalarField theta = acos(vector(1.0,0,0)&n);
    const vectorField axis = vector(1.0,0,0)^n;
    const scalarField a = cos(theta*0.5);
    const scalarField b = (axis&vector(1.0,0,0))*sin(theta*0.5);
    const scalarField c = (axis&vector(0,1.0,0))*sin(theta*0.5);
    const scalarField d = (axis&vector(0,0,1.0))*sin(theta*0.5);
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

Foam::speciesFluxConcentrationFvPatchScalarField::
speciesFluxConcentrationFvPatchScalarField
(
    const speciesFluxConcentrationFvPatchScalarField& psf
)
:
    fixedValueFvPatchScalarField(psf),
    fluxFieldName_(psf.fluxFieldName_),
    phiEName_(psf.phiEName_),
    TName_(psf.TName_),
    UName_(psf.UName_),
    chargeNumber_(psf.chargeNumber_),
    rotationMatrix_(psf.rotationMatrix_)
{}


Foam::speciesFluxConcentrationFvPatchScalarField::
speciesFluxConcentrationFvPatchScalarField
(
    const speciesFluxConcentrationFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    fluxFieldName_(psf.fluxFieldName_),
    phiEName_(psf.phiEName_),
    TName_(psf.TName_),
    UName_(psf.UName_),
    chargeNumber_(psf.chargeNumber_),
    rotationMatrix_(psf.rotationMatrix_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesFluxConcentrationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();
    
    const scalarField delta = 1.0/this->patch().deltaCoeffs();

    //const scalar t = db().time().timeOutputValue();
    //const scalar stoichCoeff = stoichCoeff_->value(t);
    //const scalar nElectrons = nElectrons_->value(t);
    
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

    // Return const access to patch values of the species flux field
    const volVectorField& speciesFlux =
         db().lookupObject<volVectorField>(fluxFieldName_+speciesName);

    const fvPatchField<vector>& speciesFluxPatch =
        patch().patchField<volVectorField, vector>(speciesFlux);
    
    // Return const access to patch values of the velocity field
    const volVectorField& U =
         db().lookupObject<volVectorField>(UName_);

    const fvPatchField<vector>& Up =
        patch().patchField<volVectorField, vector>(U);
    
    //const tmp<vectorField> tjpI =
    //      patch().patchInternalField(j.internalField());
    //const vectorField& jpI = tjpI();
    //const vectorField jMean = (jp+jpI)*0.5;

    //const tmp<vectorField>& tn = patch().nf();
    //const vectorField& n = tn();

    vectorField speciesFluxLocal = speciesFluxPatch; 
    vectorField ULocal = Up; 

    forAll(rotationMatrix_, i)
    {
        speciesFluxLocal[i] = rotationMatrix_[i]&speciesFluxPatch[i];
        ULocal[i] = rotationMatrix_[i]&Up[i];
    }
     
    const scalarField speciesFluxNormal = speciesFluxLocal&vector(1.0,0,0);
    const scalarField UNormal = ULocal&vector(1.0,0,0);
    
    if (speciesFlux.dimensions() == dimMoles/dimArea/dimTime)
    {
        const volScalarField& Deff =
            db().lookupObject<volScalarField>("Deff_" + speciesName);

        const fvPatchField<scalar>& Dp =
            patch().patchField<volScalarField, scalar>(Deff);

        const tmp<scalarField> tmup = Dp/(R*Tp);
        const scalarField& mup = tmup(); 

        const tmp<scalarField> tCpI = this->patchInternalField(); 
        const scalarField& CpI = tCpI();
        
        //const tmp<scalarField> tCp = 
        scalarField Cp = 
            //(this->patchInternalField()*Dp/delta - speciesFluxNormal)
            (CpI*Dp/delta - speciesFluxNormal)
           /(Dp/delta + double(chargeNumber_)*mup*F*(phiEp-phiEpI)/delta - UNormal);
           ///(Dp/delta - UNormal);
        //scalarField Cp = tCp();
        //forAll(Cp, faceI)
        //{
        //    if (Cp[faceI] < 0.0)
        //        Cp[faceI] = 0.0;
        //}
        
        //const tmp<scalarField> tCp = (this->patchInternalField()*Dp/delta + double(stoichCoeff)
        //        *speciesFluxNormal/(double(nElectrons)*F))
        //        /(Dp/delta + double(chargeNumber_)*mup*F*(phiEp-phiEpI)/delta);
        //const scalarField& Cp = tCp();

        operator==
        (
            (this->patchInternalField()*Dp/delta - speciesFluxNormal)
           /(Dp/delta + double(chargeNumber_)*mup*F*(phiEp-phiEpI)/delta - UNormal)
        );            
        Info << "Patch name: " << this->patch().name() << endl;
        Info << "Name of Species at boundary patch: " << speciesName << endl;
        Info << "Patch-cell distance: " << delta[0] << endl;
        Info << "Mobility: " << mup[0] << endl;
        Info << "Charge number: " << chargeNumber_ << endl;
        Info << "Diffusion coefficient: " << Dp[0] << endl;
        Info << "Patch internal concentration: " << CpI[0] << endl;
        Info << "Patch concentration: " << Cp[0] << endl;
        Info << "Patch electric potential: " << phiEp[0] << endl;
        Info << "Patch internal electric potential: " << phiEpI[0] << endl;
        Info << "Patch species flux: " << speciesFluxNormal[0] << endl;
        Info << "Patch normal velocity: " << UNormal[3] << endl;
    }
    else
    {
        FatalErrorIn("speciesFluxConcentrationFvPatchScalarField::updateCoeffs()")
            << "dimensions of species flux field are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::speciesFluxConcentrationFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>
    (
        os, 
        "fluxField", 
        "N_", 
        fluxFieldName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "electricPotentialField", 
        "phiE", 
        phiEName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "temperatureField", 
        "T", 
        TName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "velocityField", 
        "U", 
        UName_
    );
    //stoichCoeff_->writeData(os);
    //nElectrons_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//void Foam::speciesFluxConcentrationFvPatchScalarField::operator=
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
        speciesFluxConcentrationFvPatchScalarField
    );
}

// ************************************************************************* //

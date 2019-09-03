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

#include "currentConcentrationPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Foam::currentConcentrationPotentialFvPatchScalarField::
//currentConcentrationPotentialFvPatchScalarField
//(
//    const fvPatch& p,
//    const DimensionedField<scalar, volMesh>& iF
//)
//:
//    fixedValueFvPatchScalarField(p, iF),
//    currentDensityName_("i"),
//    diffCurrDensityName_("id"),
//    kappaName_("kappa"),
//    TName_("T"),
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


Foam::currentConcentrationPotentialFvPatchScalarField::
currentConcentrationPotentialFvPatchScalarField
(
    const currentConcentrationPotentialFvPatchScalarField& ccppsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ccppsf, p, iF, mapper),
    currentDensityName_(ccppsf.currentDensityName_),
    diffCurrDensityName_(ccppsf.diffCurrDensityName_),
    kappaName_(ccppsf.kappaName_),
    rotationMatrix_(ccppsf.rotationMatrix_)
{
}


//Foam::currentConcentrationPotentialFvPatchScalarField::
//currentConcentrationPotentialFvPatchScalarField
//(
//    const fvPatch& p,
//    const DimensionedField<scalar, volMesh>& iF,
//    const dictionary& dict
//)
//:
//    fixedValueFvPatchScalarField(p, iF),
//    currentDensityName_(dict.lookupOrDefault<word>("currentDensityField", "i")),
//    diffCurrDensityName_
//    (
//        dict.lookupOrDefault<word>("diffusionCurrentDensityField", "id")
//    ),
//    kappaName_(dict.lookupOrDefault<word>("conductivityField", "kappa")),
//    TName_(dict.lookupOrDefault<word>("temperatureField","T")),
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

Foam::currentConcentrationPotentialFvPatchScalarField::
currentConcentrationPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    currentDensityName_(dict.lookup("currentDensityField")),
    diffCurrDensityName_(dict.lookup("diffusionCurrentDensityField")),
    kappaName_(dict.lookup("conductivityField")),
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

Foam::currentConcentrationPotentialFvPatchScalarField::
currentConcentrationPotentialFvPatchScalarField
(
    const currentConcentrationPotentialFvPatchScalarField& ccppsf
)
:
    fixedValueFvPatchScalarField(ccppsf),
    currentDensityName_(ccppsf.currentDensityName_),
    diffCurrDensityName_(ccppsf.diffCurrDensityName_),
    kappaName_(ccppsf.kappaName_),
    rotationMatrix_(ccppsf.rotationMatrix_)
{
}

Foam::currentConcentrationPotentialFvPatchScalarField::
currentConcentrationPotentialFvPatchScalarField
(
    const currentConcentrationPotentialFvPatchScalarField& ccppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ccppsf, iF),
    currentDensityName_(ccppsf.currentDensityName_),
    diffCurrDensityName_(ccppsf.diffCurrDensityName_),
    kappaName_(ccppsf.kappaName_),
    rotationMatrix_(ccppsf.rotationMatrix_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::currentConcentrationPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalarField delta = 1.0/this->patch().deltaCoeffs();

    //const volScalarField& phiE =
    //     db().lookupObject<volScalarField>
    //     (this->dimensionedInternalField().name()); 

    //const fvPatchField<scalar>& phiEp =
    //    patch().patchField<volScalarField, scalar>(phiE);
 
    //const tmp<scalarField> tphiEpI =
    //      patch().patchInternalField(phiE.internalField());
    //const scalarField& phiEpI = tphiEpI();

    // Return const access to patch values of current density field
    const volVectorField& currentDensity =
         db().lookupObject<volVectorField>(currentDensityName_);
    const fvPatchField<vector>& currentDensityPatch =
        patch().patchField<volVectorField, vector>(currentDensity);    
    //const tmp<vectorField> tjpI =
    //      patch().patchInternalField(j.internalField());
    //const vectorField& jpI = tjpI();
    //const vectorField jMean = (jp+jpI)*0.5;

    const volVectorField& diffCurrDensity =
         db().lookupObject<volVectorField>(diffCurrDensityName_);
    const fvPatchField<vector>& diffCurrDensityPatch =
        patch().patchField<volVectorField, vector>(diffCurrDensity);    
    //const tmp<vectorField> tjdpI =
    //      patch().patchInternalField(jd.internalField());
    //const vectorField& jdpI = tjdpI();
    //const vectorField jdMean = (jdp+jdpI)*0.5;
        

     const volScalarField& kappa =
        db().lookupObject<volScalarField>(kappaName_);
    const fvPatchField<scalar>& kappap =
        patch().patchField<volScalarField, scalar>(kappa);
    //const tmp<scalarField> tkappapI =
    //      patch().patchInternalField(kappa.internalField());
    //const scalarField& kappapI = tkappapI();

    // Transform vectors into the local coordinate system

    //const tmp<vectorField>& tn = patch().nf();
    //vectorField n = tn();

    vectorField currentDensityLocal = currentDensityPatch; 
    vectorField diffCurrDensityLocal = diffCurrDensityPatch; 

    forAll(rotationMatrix_, i)
    {
        currentDensityLocal[i] = rotationMatrix_[i]&currentDensityPatch[i];
        diffCurrDensityLocal[i] = rotationMatrix_[i]&diffCurrDensityPatch[i];
    }
     
    const scalarField currentDensityNormal = currentDensityLocal&vector(1.0,0,0);
    const scalarField diffCurrDensityNormal = diffCurrDensityLocal&vector(1.0,0,0);
    //const scalarField kappaMean = (kappap+kappapI)*0.5;
    //const scalarField kappaMean = kappap*kappapI*2.0/(kappap+kappapI);
    //const scalarField kappaMean = sqrt(kappap*kappapI);
    //tmp<scalarField> ttest = this->patchInternalField() 
    //                + (diffCurrDensityNormal - currentDensityNormal)*delta/kappa;
    //scalarField test = ttest();

    if (currentDensity.dimensions() == dimCurrent/dimArea)
    {
        if (min(kappap) > 0.0)
        {
            scalarField::operator=
            (
                this->patchInternalField() 
                    + (diffCurrDensityNormal - currentDensityNormal)*delta/kappap
            );
            //Info << "Patch name: " << this->patch().name() << endl;
            //Info << "Patch conductivity: " << kappap[0] << endl;
            //Info << "Patch-cell distance: " << delta[0] << endl;
            //Info << "Patch normal current density (local coordinates): " << currentDensityNormal[0] << endl;
            //Info << "Patch normal diffusion current density (local coordinates): " << diffCurrDensityNormal[0] << endl;
            //Info << "Patch concentration value: " << test[0] << endl;
            //Info << "Patch axis vector: " << axis[0] << endl;
            //Info << "Patch diffusion current density (local coordinates): " << jdn[0] << endl;
            //Info << "Patch normal: " << n[0] << endl; 
        }
        else
        {
            scalarField::operator=
            (
                this->patchInternalField() 
            );
        }
    }
    else
    {
        FatalErrorIn("currentConcentrationPotentialFvPatchScalarField::updateCoeffs()")
            << "dimensions of current density field are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::currentConcentrationPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>
    (
        os, "currentDensityField", "currentDensityField", currentDensityName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "diffusionCurrentDensityField", 
        "diffusionCurrentDensityField", 
        diffCurrDensityName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "conductivityField", 
        "conductivityField", 
        kappaName_
    );
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
void Foam::currentConcentrationPotentialFvPatchScalarField::operator=
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
        currentConcentrationPotentialFvPatchScalarField
    );
}

// ************************************************************************* //

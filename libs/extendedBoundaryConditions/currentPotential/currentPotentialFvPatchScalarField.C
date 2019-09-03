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

#include "currentPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::currentPotentialFvPatchScalarField::
currentPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    iName_("i"),
    idName_("id"),
    rotationMatrix_(patch().nf()->size())
{
    // Set up rotation matrix
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


Foam::currentPotentialFvPatchScalarField::
currentPotentialFvPatchScalarField
(
    const currentPotentialFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    iName_(psf.iName_),
    idName_(psf.idName_),
    rotationMatrix_(psf.rotationMatrix_)
{}


Foam::currentPotentialFvPatchScalarField::
currentPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    iName_(dict.lookupOrDefault<word>("currentDensityField", "i")),
    idName_(dict.lookupOrDefault<word>("diffusionCurrentDensityField", "id")),
    rotationMatrix_(patch().nf()->size())
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    // Set up rotation matrix
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


Foam::currentPotentialFvPatchScalarField::
currentPotentialFvPatchScalarField
(
    const currentPotentialFvPatchScalarField& psf
)
:
    fixedValueFvPatchScalarField(psf),
    iName_(psf.iName_),
    idName_(psf.idName_),
    rotationMatrix_(psf.rotationMatrix_)
{
}

Foam::currentPotentialFvPatchScalarField::
currentPotentialFvPatchScalarField
(
    const currentPotentialFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    iName_(psf.iName_),
    idName_(psf.idName_),
    rotationMatrix_(psf.rotationMatrix_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::currentPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
   
    // Physical constants
    //const scalar& R = constant::physicoChemical::R.value();
    //const scalar& F = constant::physicoChemical::F.value();
    
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
    const volVectorField& i =
         db().lookupObject<volVectorField>(iName_);
    const fvPatchField<vector>& ip =
        patch().patchField<volVectorField, vector>(i);    

    const volVectorField& id =
         db().lookupObject<volVectorField>(idName_);
    const fvPatchField<vector>& idp =
        patch().patchField<volVectorField, vector>(id);    
        
     const volScalarField& kappa =
        db().lookupObject<volScalarField>("kappa");
    const fvPatchField<scalar>& kappap =
        patch().patchField<volScalarField, scalar>(kappa);
    const tmp<scalarField> tkappapI =
          patch().patchInternalField(kappa.internalField());
    const scalarField& kappapI = tkappapI();

    // Transform vectors into the local coordinate system
    vectorField iLocal = ip; 
    vectorField idLocal = idp; 
    forAll(rotationMatrix_, i)
    {
        
            iLocal[i] = rotationMatrix_[i]&ip[i];
            idLocal[i] = rotationMatrix_[i]&idp[i];
    }
     
    const scalarField iNormal = iLocal&vector(1.0,0,0);
    const scalarField idNormal = idLocal&vector(1.0,0,0);
    //const scalarField kappaMean = (kappap+kappapI)/2.0;
    //const scalarField kappaMean = kappap*kappapI*2.0/(kappap+kappapI);
    //const scalarField kappaMean = sqrt(kappap*kappapI);

    //const Field<scalar>& magS = patch().magSf();
              
    if (i.dimensions() == dimCurrent/dimArea)
    {
        if (min(kappap) > 0.0)
        {
            scalarField::operator=
            (
                this->patchInternalField() + (idNormal - iNormal)*delta/kappap
            );
            Info << "Patch name: " << this->patch().name() << endl;
            Info << "Patch conductivity: " << kappap[0] << endl;
            Info << "Patch-cell distance: " << delta[0] << endl;
            Info << "Patch current density (local coordinates): " << iNormal[0] << endl;
            Info << "Patch diffusion current density (local coordinates): " << idNormal[0] << endl;
        }
        else
        {
            FatalErrorIn("currentPotentialFvPatchScalarField::updateCoeffs()")
                << "Ionic conductivity kappa must be > 0"
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("currentPotentialFvPatchScalarField::updateCoeffs()")
            << "dimensions of j are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::currentPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "currentDensityField", "i", iName_);
    writeEntryIfDifferent<word>(os, "diffusionCurrentDensityField", "id", idName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
void Foam::currentPotentialFvPatchScalarField::operator=
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
        currentPotentialFvPatchScalarField
    );
}

// ************************************************************************* //

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

#include "thixoShearStressSlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thixoShearStressSlipVelocityFvPatchVectorField::
thixoShearStressSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    modelName_("Navier"),
    nuName_("nu"),
    lambdaName_("lambda"),
    alpha_(1.0),
    beta_(1.0),
    gamma_(1.0),
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


Foam::thixoShearStressSlipVelocityFvPatchVectorField::
thixoShearStressSlipVelocityFvPatchVectorField
(
    const thixoShearStressSlipVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    modelName_(pvf.modelName_),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    beta_(pvf.beta_),
    gamma_(pvf.gamma_),
    rho_(pvf.rho_)
{}


Foam::thixoShearStressSlipVelocityFvPatchVectorField::
thixoShearStressSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    modelName_(dict.lookupOrDefault<word>("modelName","Navier")),
    nuName_(dict.lookupOrDefault<word>("nuName","nu")),
    lambdaName_(dict.lookupOrDefault<word>("lambdaName","lambda")),
    alpha_(readScalar(dict.lookup("alpha"))),
    beta_(readScalar(dict.lookup("beta"))),
    gamma_(readScalar(dict.lookup("gamma"))),
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


Foam::thixoShearStressSlipVelocityFvPatchVectorField::
thixoShearStressSlipVelocityFvPatchVectorField
(
    const thixoShearStressSlipVelocityFvPatchVectorField& pvf 
)
:
    fixedValueFvPatchVectorField(pvf),
    modelName_(pvf.modelName_),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    beta_(pvf.beta_),
    gamma_(pvf.gamma_),
    rho_(pvf.rho_)
{}


Foam::thixoShearStressSlipVelocityFvPatchVectorField::
thixoShearStressSlipVelocityFvPatchVectorField
(
    const thixoShearStressSlipVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    modelName_(pvf.modelName_),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    beta_(pvf.beta_),
    gamma_(pvf.gamma_),
    rho_(pvf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thixoShearStressSlipVelocityFvPatchVectorField::updateCoeffs()
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
    const volScalarField& lambda = 
        db().lookupObject<volScalarField>(lambdaName_); 
    const scalarField& lambdap =
        patch().patchField<volScalarField, scalar>(lambda);
    const volScalarField& nuPrev = nu.prevIter();
    const scalarField& nuPrevp =
        patch().patchField<volScalarField, scalar>(nuPrev);
    
    tmp<vectorField> tn = patch().nf();
    const vectorField& n = tn();

    const vectorField tUprevp = Uprevp - n*(n & Uprevp);
    const vectorField tUprevpI = UprevpI() - n*(n & UprevpI());
    //const vectorField nUgradPrev = (tUprevp - tUprevpI)*patch().deltaCoeffs();

    const vectorField g = (tUprevp - tUprevpI)*patch().deltaCoeffs();
    
    //const scalarField d = alpha_ * pow(gamma_,lambdap)
    scalarField d = nuPrevp*0.0;
    if(modelName_ == "Navier")
    {
        d = max
            (
                alpha_ * pow(max(mag(g),VSMALL), (beta_-1.0))
                * pow(nuPrevp*rho_.value(), beta_)
                * patch().deltaCoeffs(),
                0.0
            );
    }
    else if(modelName_ == "Hatzikiriakos")
    {
        d = max
            (
                alpha_ * sinh(beta_*mag(nuPrevp*g)) *  nuPrevp
                / max(mag(nuPrevp*g),VSMALL) * patch().deltaCoeffs(),
                0.0
            );
    }
    else if(modelName_ == "asymptotic")
    {
        d = max
            (
                alpha_ * log(1 + beta_*mag(nuPrevp*g)) *  nuPrevp
                / max(mag(nuPrevp*g),VSMALL) * patch().deltaCoeffs(),
                0.0
            );
    }
    else
    {
        FatalErrorIn
        (
            "thixoShearStressSlipVelocityFvPatchVectorField::updateCoeffs()"
        )   << "Unknown slip velocity model "
            << modelName_ << endl << endl
            << "Valid models are : " << endl
            << "Navier" << endl
            << "Hatzikiriakos" << endl
            << "asymptotic" << exit(FatalError);
    }

    d *= exp(-gamma_ * lambdap);
    vectorField Up = d/(d+1.0)*tUprevpI;

    operator==(Up);


    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::thixoShearStressSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>
    (
        os, 
        "modelName", 
        "Navier", 
        modelName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "nuName", 
        "nu", 
        nuName_
    );
    writeEntryIfDifferent<word>
    (
        os, 
        "lambdaName", 
        "lambda", 
        lambdaName_
    );
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//void Foam::thixoShearStressSlipVelocityFvPatchVectorField::operator=
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
        thixoShearStressSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

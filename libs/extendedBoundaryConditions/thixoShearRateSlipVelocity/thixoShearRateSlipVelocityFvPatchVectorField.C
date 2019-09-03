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

#include "thixoShearRateSlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thixoShearRateSlipVelocityFvPatchVectorField::
thixoShearRateSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    nuName_("nu"),
    lambdaName_("lambda"),
    alpha_(1.0),
    gamma_(1.0),
    beta_(1.0),
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


Foam::thixoShearRateSlipVelocityFvPatchVectorField::
thixoShearRateSlipVelocityFvPatchVectorField
(
    const thixoShearRateSlipVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    gamma_(pvf.gamma_),
    beta_(pvf.beta_),
    rho_(pvf.rho_)
{}


Foam::thixoShearRateSlipVelocityFvPatchVectorField::
thixoShearRateSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    nuName_(dict.lookupOrDefault<word>("nuName","nu")),
    lambdaName_(dict.lookupOrDefault<word>("lambdaName","lambda")),
    alpha_(readScalar(dict.lookup("alpha"))),
    gamma_(readScalar(dict.lookup("gamma"))),
    beta_(readScalar(dict.lookup("beta"))),
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


Foam::thixoShearRateSlipVelocityFvPatchVectorField::
thixoShearRateSlipVelocityFvPatchVectorField
(
    const thixoShearRateSlipVelocityFvPatchVectorField& pvf 
)
:
    fixedValueFvPatchVectorField(pvf),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    gamma_(pvf.gamma_),
    beta_(pvf.beta_),
    rho_(pvf.rho_)
{}


Foam::thixoShearRateSlipVelocityFvPatchVectorField::
thixoShearRateSlipVelocityFvPatchVectorField
(
    const thixoShearRateSlipVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    gamma_(pvf.gamma_),
    beta_(pvf.beta_),
    rho_(pvf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thixoShearRateSlipVelocityFvPatchVectorField::updateCoeffs()
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
    const vectorField nUgradPrev = (tUprevp - tUprevpI)*patch().deltaCoeffs();
    
    const scalarField d = alpha_*exp(-gamma_*lambdap)
        * pow(max(mag(nUgradPrev),VSMALL), (beta_-1.0))
        * patch().deltaCoeffs();
    
    vectorField Up = d/(d+1.0)*tUprevpI;

    operator==(Up);


    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::thixoShearRateSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
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
    os.writeKeyword("alpha") << alpha_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ 
        << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//void Foam::thixoShearRateSlipVelocityFvPatchVectorField::operator=
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
        thixoShearRateSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

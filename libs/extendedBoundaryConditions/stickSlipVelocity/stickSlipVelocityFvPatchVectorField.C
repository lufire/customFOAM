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

#include "stickSlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stickSlipVelocityFvPatchVectorField::
stickSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    nuName_("nu"),
    lambdaName_("lambda"),
    alpha_(1.0),
    beta_(0.0),
    gamma_(0.0),
    sg_(0.0),
    rho_(0.0),
    tau0_(0.0),
    Ustar_(0.0)
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
    const dimensionedScalar rho = transDict.lookup("rho");
    rho_ = rho.value();
    const word transModel = transDict.lookup("transportModel");
    const word transModelCoeffs = transModel + "Coeffs";
    const dictionary& transModelDict = transDict.subDict(transModelCoeffs);
    const dimensionedScalar tau0 = transModelDict.lookup("tau0");
    tau0_ = tau0.value();
}


Foam::stickSlipVelocityFvPatchVectorField::
stickSlipVelocityFvPatchVectorField
(
    const stickSlipVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    beta_(pvf.beta_),
    gamma_(pvf.gamma_),
    sg_(pvf.sg_),
    rho_(pvf.rho_),
    tau0_(pvf.tau0_),
    Ustar_(pvf.Ustar_)
{}


Foam::stickSlipVelocityFvPatchVectorField::
stickSlipVelocityFvPatchVectorField
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
    beta_(readScalar(dict.lookup("beta"))),
    gamma_(readScalar(dict.lookup("gamma"))),
    sg_(readScalar(dict.lookup("sg"))),
    rho_(dict.lookupOrDefault<scalar>("rho", 0.0)),
    tau0_(dict.lookupOrDefault<scalar>("tau0", 0.0)),
    //rho_
    //(
    //    dict.lookupOrDefault<dimensionedScalar>
    //    (
    //        "rho",
    //        dimensionedScalar("zero", dimMass/dimVol, 0.0)
    //    )
    //),
    //tau0_
    //(
    //    dict.lookupOrDefault<dimensionedScalar>
    //    (
    //        "tau0",
    //        dimensionedScalar("zero", dimViscosity/dimTime, 0.0)
    //    )
    //),
    Ustar_(readScalar(dict.lookup("Ustar")))
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
    if(rho_ == 0.0)
    {
        const dimensionedScalar rho = transDict.lookup("rho");
        rho_ = rho.value();
    }
    if(tau0_ == 0.0)
    {
        const word transModel = transDict.lookup("transportModel");
        const word transModelCoeffs = transModel + "Coeffs";
        const dictionary& transModelDict = transDict.subDict(transModelCoeffs);
        const dimensionedScalar tau0 = transModelDict.lookup("tau0");
        tau0_ = tau0.value();
    }
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::stickSlipVelocityFvPatchVectorField::
stickSlipVelocityFvPatchVectorField
(
    const stickSlipVelocityFvPatchVectorField& pvf 
)
:
    fixedValueFvPatchVectorField(pvf),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    beta_(pvf.beta_),
    gamma_(pvf.gamma_),
    sg_(pvf.sg_),
    rho_(pvf.rho_),
    tau0_(pvf.tau0_),
    Ustar_(pvf.Ustar_)
{}


Foam::stickSlipVelocityFvPatchVectorField::
stickSlipVelocityFvPatchVectorField
(
    const stickSlipVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    nuName_(pvf.nuName_),
    lambdaName_(pvf.lambdaName_),
    alpha_(pvf.alpha_),
    beta_(pvf.beta_),
    gamma_(pvf.gamma_),
    sg_(pvf.sg_),
    rho_(pvf.rho_),
    tau0_(pvf.tau0_),
    Ustar_(pvf.Ustar_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stickSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const volVectorField& U =
         db().lookupObject<volVectorField>("U");
    const volVectorField& Uprev = U.prevIter();
    const tmp<vectorField> UpI =
          patch().patchInternalField(U.internalField());
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
    const scalarField& nup =
        patch().patchField<volScalarField, scalar>(nu);
    
    tmp<vectorField> tn = patch().nf();
    const vectorField& n = tn();

    const vectorField tUprevp = Uprevp - n*(n & Uprevp);
    const vectorField tUprevpI = UprevpI() - n*(n & UprevpI());
    const vectorField tUpI = UpI() - n*(n & UpI());
    const vectorField nUgradPrev = 
        (tUprevp - tUprevpI)*patch().deltaCoeffs();

    //const vectorField tauW = 
    //    (tUprevp - tUprevpI)*patch().deltaCoeffs()*rho_*nup;
    const scalarField magTauW = mag(nUgradPrev*rho_*nup);

    const scalarField magTUpI = mag(tUpI);
    const scalar rAlpha = 1.0/alpha_;

    const scalarField d = alpha_* pow(max(mag(nUgradPrev),VSMALL), (beta_-1.0))
        * patch().deltaCoeffs();
    const scalarField factor = d/(d + 1.0);
    vectorField Up = tUpI*0.0;
    forAll(Up, faceI)
    {
        if(magTauW[faceI] <= sg_)
        {
            Up[faceI] *= 0.0; 
        }
        else
        {
            Up[faceI] = factor[faceI]*tUprevpI[faceI];
        }
    }
    Up *= exp(-gamma_ * lambdap);

    operator==(Up);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::stickSlipVelocityFvPatchVectorField::write(Ostream& os) const
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
    os.writeKeyword("beta") << beta_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("sg") << sg_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("tau0") << tau0_ 
        << token::END_STATEMENT << nl;
    os.writeKeyword("Ustar") << Ustar_ 
        << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//void Foam::stickSlipVelocityFvPatchVectorField::operator=
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
        stickSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

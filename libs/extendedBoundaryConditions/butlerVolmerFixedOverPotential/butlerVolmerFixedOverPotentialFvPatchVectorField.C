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

#include "butlerVolmerFixedOverPotentialFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::
butlerVolmerFixedOverPotentialFvPatchVectorField
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
    iEx_(0.0),
    eqPotential_(0.0),
    electrodePotential_(0.0),
    alphaA_(1.0),
    alphaC_(1.0),
    nElectrons_(1)
{}


Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::
butlerVolmerFixedOverPotentialFvPatchVectorField
(
    const butlerVolmerFixedOverPotentialFvPatchVectorField& ptf,
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
    iEx_(ptf.iEx_),
    eqPotential_(ptf.eqPotential_),
    electrodePotential_(ptf.electrodePotential_),
    alphaA_(ptf.alphaA_),
    alphaC_(ptf.alphaC_),
    nElectrons_(ptf.nElectrons_)
{}


Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::
butlerVolmerFixedOverPotentialFvPatchVectorField
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
    //oxidantName_(dict.lookup("oxidant")),
    //reductantName_(dict.lookup("reductant")),
    //CRefOxidant_(readScalar(dict.lookup("CRefOxidant"))),
    //CRefReductant_(readScalar(dict.lookup("CRefReductant"))),
    //gammaOxidant_(readScalar(dict.lookup("gammaOxidant"))),
    //gammaReductant_(readScalar(dict.lookup("gammaReductant"))),
    //iEx_(readScalar(dict.lookup("exchangeCurrentDensity"))),
    //eqPotential_(readScalar(dict.lookup("equilibriumPotential"))),
    //electrodePotential_(readScalar(dict.lookup("electrodePotential"))),
    //alphaA_(readScalar(dict.lookup("alphaA"))),
    //alphaC_(readScalar(dict.lookup("alphaC"))),
    //nElectrons_(readLabel(dict.lookup("electrons")))
    oxidantName_("oxidant"),
    reductantName_("reductant"),
    CRefOxidant_(1.0),
    CRefReductant_(1.0),
    gammaOxidant_(1.0),
    gammaReductant_(1.0),
    iEx_(0.0),
    eqPotential_(0.0),
    electrodePotential_(0.0),
    alphaA_(1.0),
    alphaC_(1.0),
    nElectrons_(1)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::
butlerVolmerFixedOverPotentialFvPatchVectorField
(
    const butlerVolmerFixedOverPotentialFvPatchVectorField& bvcpvf
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
    iEx_(bvcpvf.iEx_),
    eqPotential_(bvcpvf.eqPotential_),
    electrodePotential_(bvcpvf.electrodePotential_),
    alphaA_(bvcpvf.alphaA_),
    alphaC_(bvcpvf.alphaC_),
    nElectrons_(bvcpvf.nElectrons_)
{}


Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::
butlerVolmerFixedOverPotentialFvPatchVectorField
(
    const butlerVolmerFixedOverPotentialFvPatchVectorField& bvcpvf,
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
    iEx_(bvcpvf.iEx_),
    eqPotential_(bvcpvf.eqPotential_),
    electrodePotential_(bvcpvf.electrodePotential_),
    alphaA_(bvcpvf.alphaA_),
    alphaC_(bvcpvf.alphaC_),
    nElectrons_(bvcpvf.nElectrons_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Const reference to temperature 
    const volScalarField& T =
        db().lookupObject<volScalarField>(TName_);
    const fvPatchField<scalar>& Tp =
        patch().patchField<volScalarField, scalar>(T);
    if(gMin(Tp) <= 0.0)
    {
        FatalErrorInFunction
            << "Temperature must be > 0.0 " << endl
            << "on patch " << this->patch().name() << "."
            << exit(FatalError);
    }
    
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();

    // Assign variables
    const dictionary& electrochemDict = 
        db().lookupObject<IOdictionary>("electrochemicalProperties");
    const dictionary electrodeDict = electrochemDict.subDict(electrodeName_);
    oxidantName_ = electrodeDict.lookupType<word>("oxidant");
    reductantName_ = electrodeDict.lookupType<word>("reductant"); 
    const word oxName = oxidantName_;
    const word redName = reductantName_;
    CRefOxidant_ = 
        electrodeDict.subDict("CRefs").lookupType<scalar>(oxName);
    CRefReductant_ = 
        electrodeDict.subDict("CRefs").lookupType<scalar>(redName);
    if (CRefOxidant_ <= 0.0 || CRefReductant_ <= 0.0)
    {
        FatalErrorInFunction
            << "Reference concentrations must be greater than zero"
            << exit(FatalError);
    }
    gammaOxidant_ = 
        electrodeDict.subDict("reactionOrders").lookupType<scalar>(oxName);
    gammaReductant_ =
        electrodeDict.subDict("reactionOrders").lookupType<scalar>(redName);
    iEx_ = electrodeDict.lookupType<scalar>("exchangeCurrentDensity");
    eqPotential_ = electrodeDict.lookupType<scalar>("equilibriumPotential");
    electrodePotential_ = 
        electrodeDict.lookupType<scalar>("potential");
    alphaA_ = electrodeDict.lookupType<scalar>("alphaA");
    alphaC_ = electrodeDict.lookupType<scalar>("alphaC");
    nElectrons_ = electrodeDict.lookupType<label>("electrons");

    // Const reference to electric potential 
    const volScalarField& phiE =
         db().lookupObject<volScalarField>(phiEName_); 
    const fvPatchField<scalar>& phiEls =
        patch().patchField<volScalarField, scalar>(phiE);
 
    // Const reference to oxidant concentration
    const volScalarField& COxidant =
        db().lookupObject<volScalarField>("C_"+oxidantName_);
    const fvPatchField<scalar>& COxp =
        patch().patchField<volScalarField, scalar>(COxidant);

    // Const reference to reductant concentration
    const volScalarField& CReductant =
        db().lookupObject<volScalarField>("C_"+reductantName_);
    const fvPatchField<scalar>& CRedp =
        patch().patchField<volScalarField, scalar>(CReductant);

    // Patch normal vectors
    const tmp<vectorField>& tn = patch().nf();
    const vectorField& n = tn();
    const scalarField& magSf = patch().magSf();
    const scalar totalMagSf = gSum(magSf);

    if (phiE.dimensions() == dimMass*dimArea/pow3(dimTime)/dimCurrent)
    {
        vectorField ip = *this * 0.0;
        scalarField redFactor = phiEls * 0.0;
        scalarField oxFactor = redFactor;
        const scalar minFactor = SMALL;
        forAll(redFactor, faceI)
        {
            redFactor[faceI] = 
                pow(CRedp[faceI]/CRefReductant_, gammaReductant_);
            oxFactor[faceI] = 
                pow(COxp[faceI]/CRefOxidant_, gammaOxidant_);
            if(redFactor[faceI] <= 0.0) 
            {
                redFactor[faceI] = minFactor;
            }
            if(oxFactor[faceI] <= 0.0)
            {
                oxFactor[faceI] = minFactor;
            }
        }
        scalarField K = R*Tp/(nElectrons_*F);
        // Initial guess of electrode potential via Tafel equation
        scalarField eta = phiEls * 0.0;
        //const scalar alphaC = 1.0-alphaA_;
        scalar iAvg = 0.0;
        scalar maxCurrent = 1e4;
        forAll(n, faceI)
        {
            //eta[faceI] = electrodePotential_ - phiEls[faceI] - eqPotential_
            //    + K[faceI]*log(oxFactor[faceI]/redFactor[faceI]);
            //Info << "faceI: " << faceI 
            //    << ", eta: " << eta[faceI] 
            //    << ", electrodePotential: " << electrodePotential_
            //    << ", phiEls: " << phiEls[faceI] 
            //    << ", oxFactor: " << oxFactor[faceI] 
            //    << ", redFactor: " << redFactor[faceI] << endl;
            //scalar i0 = iEx_
            //    * pow(oxFactor[faceI], alphaC_)
            //    * pow(redFactor[faceI], alphaA_);
            //scalar iT = i0 
            //    * (exp(alphaA_/K[faceI]*eta[faceI]) 
            //       - exp(-alphaC_/K[faceI]*eta[faceI]));
            eta[faceI] = electrodePotential_ - phiEls[faceI] - eqPotential_;
            Info << "faceI: " << faceI 
                << ", eta: " << eta[faceI] 
                << ", electrodePotential: " << electrodePotential_
                << ", phiEls: " << phiEls[faceI] 
                << ", COxp: " << COxp[faceI] 
                << ", CRedp: " << CRedp[faceI]
                << ", oxFactor: " << oxFactor[faceI] 
                << ", redFactor: " << redFactor[faceI] << endl;
            scalar iT = iEx_ 
               * (redFactor[faceI]*exp((alphaA_)/K[faceI]*eta[faceI])
                  - oxFactor[faceI]*exp(-alphaC_)/K[faceI]*eta[faceI]);
            iT = (mag(iT)<maxCurrent)?iT:Foam::sign(iT)*maxCurrent;
            ip[faceI] = n[faceI]*-1.0*iT;
            Info << "i["<<faceI<<"] = " << ip[faceI] << endl;
            iAvg += magSf[faceI]*iT;
        }
        iAvg /= totalMagSf;
        Info << "iAvg: " << iAvg << endl;
        operator==(ip);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiEName_ << " are not correct" << endl
            << "on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "electrodeName", electrodeName_);
    writeEntryIfDifferent<word>(os, "phiE", "phiE", phiEName_);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::butlerVolmerFixedOverPotentialFvPatchVectorField::operator=
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
        butlerVolmerFixedOverPotentialFvPatchVectorField
    );
}

// ************************************************************************* //

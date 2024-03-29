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

#include "butlerVolmerCurrentFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::butlerVolmerCurrentFvPatchVectorField::
butlerVolmerCurrentFvPatchVectorField
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


Foam::butlerVolmerCurrentFvPatchVectorField::
butlerVolmerCurrentFvPatchVectorField
(
    const butlerVolmerCurrentFvPatchVectorField& ptf,
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


Foam::butlerVolmerCurrentFvPatchVectorField::
butlerVolmerCurrentFvPatchVectorField
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
    gammaOxidant_ = 
        electrodeDict.subDict("reactionOrders").lookupType<scalar>(oxName);
    gammaReductant_ =
        electrodeDict.subDict("reactionOrders").lookupType<scalar>(redName);
    iEx_ = electrodeDict.lookupType<scalar>("exchangeCurrentDensity");
    eqPotential_ = electrodeDict.lookupType<scalar>("equilibriumPotential");
    electrodePotential_ = 
        electrodeDict.lookupType<scalar>("electrodePotential");
    alphaA_ = electrodeDict.lookupType<scalar>("alphaA");
    alphaC_ = electrodeDict.lookupType<scalar>("alphaC");
    nElectrons_ = electrodeDict.lookupType<label>("electrons");
}


Foam::butlerVolmerCurrentFvPatchVectorField::
butlerVolmerCurrentFvPatchVectorField
(
    const butlerVolmerCurrentFvPatchVectorField& bvcpvf
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


Foam::butlerVolmerCurrentFvPatchVectorField::
butlerVolmerCurrentFvPatchVectorField
(
    const butlerVolmerCurrentFvPatchVectorField& bvcpvf,
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

void Foam::butlerVolmerCurrentFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();

    // Return const access to internal cell values adjacent to
    // patch of electric field
    const volScalarField& phiE =
         db().lookupObject<volScalarField>(phiEName_); 
    const fvPatchField<scalar>& phiEp =
        patch().patchField<volScalarField, scalar>(phiE);
 
    // Return const access to internal cell values adjacent to 
    // patch of temperature
    const volScalarField& T =
        db().lookupObject<volScalarField>(TName_);
    const fvPatchField<scalar>& Tp =
        patch().patchField<volScalarField, scalar>(T);
         
    // Return const access to internal cell values adjacent to patch of oxidant
    // mass fraction
    const volScalarField& COxidant =
        db().lookupObject<volScalarField>("C_"+oxidantName_);
    const fvPatchField<scalar>& COxidantp =
        patch().patchField<volScalarField, scalar>(COxidant);

    // Return const access to internal cell values adjacent to patch of 
    // reductant mass fraction    
    const volScalarField& CReductant =
        db().lookupObject<volScalarField>("C_"+reductantName_);
    const fvPatchField<scalar>& CReductantp =
        patch().patchField<volScalarField, scalar>(CReductant);
 
    const tmp<vectorField>& tn = patch().nf();
    const vectorField& n = tn();
//    const Field<scalar>& magS = patch().magSf();
           
    vectorField ip(phiEp.size(), vector::zero);
        
    scalar maxCurrent = 1e4;
    if (phiE.dimensions() == dimMass*dimArea/pow3(dimTime)/dimCurrent)
    {
        forAll(ip, faceI)
        {
            scalar redFactor = 
                pow(CReductantp[faceI]/CRefReductant_, gammaReductant_);
            scalar oxFactor = 
                pow(COxidantp[faceI]/CRefOxidant_,gammaOxidant_);
            scalar K = R*Tp[faceI]/(nElectrons_*F);
            scalar i0 = iEx_*pow(redFactor,alphaC_)*pow(oxFactor,alphaA_);
            scalar overPotential = 
                electrodePotential_ - phiEp[faceI] - eqPotential_;
            //if(mag(overPotential) > mag(eqPotential_))
            //{
            //    overPotential = Foam::sign(overPotential)*mag(eqPotential_);
            //}
            //scalar iT = iEx_
            //   *(redFactor*Foam::exp(alphaA_/K*overPotential)
            //     -oxFactor*Foam::exp(-alphaC_/K*overPotential));
            scalar iT = i0
               *(exp(alphaA_/K
                     *(overPotential + K*log(redFactor/oxFactor)))
                 -exp(-(1-alphaA_)/K
                      *(overPotential + K*log(redFactor/oxFactor))));

                 //  Foam::exp(alphaA_*F/R/Tp[faceI]
                 // *(electrodePotential_-phiEp[faceI]-eqPotential_))                  
                 // -Foam::exp(-alphaC_*F/R/Tp[faceI]
                 // *(electrodePotential_-phiEp[faceI]-eqPotential_))


            //ip[faceI] =  n[faceI]*pow(YReductantp[faceI]/CRefReductant_,gammaReductant_);
            //iT = (mag(iT)<maxCurrent)?iT:Foam::sign(iT)*maxCurrent;
            iT = (mag(iT)<maxCurrent)?iT:Foam::sign(iT)*maxCurrent;
            ip[faceI] = iT*n[faceI]*-1.0;

            Info << "faceI" << faceI
            << ", CReductantp = " << CReductantp[faceI]                       
            << ", COxidantp = " << COxidantp[faceI]
            << ", CRefReductant = " << CRefReductant_
            << ", CRefOxidant = " << CRefOxidant_
            << ", phiE = " << phiEp[faceI] 
            << ", overPotential = " << overPotential 
            << ", oxFactor = " << oxFactor 
            << ", redFactor = " << redFactor
            << "ip = " << ip[faceI] << endl; 
        }
        operator==(ip);
//             
    }
    else
    {
        FatalErrorIn("butlerVolmerCurrentFvPatchVectorField::updateCoeffs()")
            << "dimensions of phiE are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
//    }
}

void Foam::butlerVolmerCurrentFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "phiE", phiEName_);
    writeEntry(os, "T", TName_);
    writeEntry(os, "oxidant", oxidantName_);
    writeEntry(os, "reductant", reductantName_);
    writeEntry(os, "CRefOxidant", CRefOxidant_);
    writeEntry(os, "CRefReductant", CRefReductant_);
    writeEntry(os, "gammaOxidant", gammaOxidant_);
    writeEntry(os, "gammaReductant", gammaReductant_);
    writeEntry(os, "exchangeCurrentDensity", iEx_);
    writeEntry(os, "equilibriumPotential", eqPotential_);
    writeEntry(os, "electrodePotential", electrodePotential_);
    writeEntry(os, "alphaA", alphaA_);
    writeEntry(os, "alphaC", alphaC_);
    writeEntry(os, "electrons", nElectrons_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::butlerVolmerCurrentFvPatchVectorField::operator=
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
        butlerVolmerCurrentFvPatchVectorField
    );
}

// ************************************************************************* //

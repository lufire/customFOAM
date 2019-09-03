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

#include "butlerVolmerPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::butlerVolmerPotentialFvPatchScalarField::
butlerVolmerPotentialFvPatchScalarField
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
    alphaA_(1.0),
    alphaC_(1.0)
{}


Foam::butlerVolmerPotentialFvPatchScalarField::
butlerVolmerPotentialFvPatchScalarField
(
    const butlerVolmerPotentialFvPatchScalarField& ptf,
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
    alphaA_(ptf.alphaA_),
    alphaC_(ptf.alphaC_)
{}


Foam::butlerVolmerPotentialFvPatchScalarField::
butlerVolmerPotentialFvPatchScalarField
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
    alphaA_(dict.lookupOrDefault<scalar>("alphaA", 1.0)),
    alphaC_(dict.lookupOrDefault<scalar>("alphaC", 1.0))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::butlerVolmerPotentialFvPatchScalarField::
butlerVolmerPotentialFvPatchScalarField
(
    const butlerVolmerPotentialFvPatchScalarField& bvppsf
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
    alphaA_(bvppsf.alphaA_),
    alphaC_(bvppsf.alphaC_)
{}


Foam::butlerVolmerPotentialFvPatchScalarField::
butlerVolmerPotentialFvPatchScalarField
(
    const butlerVolmerPotentialFvPatchScalarField& bvppsf,
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
    alphaA_(bvppsf.alphaA_),
    alphaC_(bvppsf.alphaC_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::butlerVolmerPotentialFvPatchScalarField::bvNewtonMethodFunction
(
    const scalar& overPotential,
    const label& faceI
)
{

}

void Foam::butlerVolmerPotentialFvPatchScalarField::bvNewtonMethodDerivative
(
    const scalar& overPotential
)
{

}
*/
    

void Foam::butlerVolmerPotentialFvPatchScalarField::updateCoeffs()
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
         db().lookupObject<volScalarField>
         (this->dimensionedInternalField().name()); 

    const fvPatchField<scalar>& phiEp =
        patch().patchField<volScalarField, scalar>(phiE);
 
    //const tmp<scalarField> tphiEpI =
    //      patch().patchInternalField(phiE.internalField());
          

    //const scalarField& phiEpI = tphiEpI();


    // Return const access to patch values of current density field
    const volVectorField& j =
         db().lookupObject<volVectorField>(jName_);

    const fvPatchField<vector>& jp =
        patch().patchField<volVectorField, vector>(j);
 
    //const tmp<vectorField> tjpI =
    //      patch().patchInternalField(j.internalField());
          
    //const vectorField& jpI = tjpI();
    
    // Return const access to internal cell values adjacent to 
    // patch of temperature
    const volScalarField& T =
        db().lookupObject<volScalarField>(TName_);
        
    const fvPatchField<scalar>& Tp =
        patch().patchField<volScalarField, scalar>(T);
 
    //const tmp<scalarField> tTpI =
    //    patch().patchInternalField(T.internalField());

    //const scalarField& TpI = tTpI();

    //const dictionary& electrochemicalProperties =
    //    db().lookupObject<IOdictionary>("electrochemicalProperties");
    //electrochemicalProperties.subDict(electrodeName_).lookup("");

    // Return const access to internal cell values adjacent to patch of oxidant
    // mass fraction
    const volScalarField& COxidant =
        db().lookupObject<volScalarField>(electrodeName_+"Oxidant");
               

    const fvPatchField<scalar>& COxidantp =
        patch().patchField<volScalarField, scalar>(COxidant);

//    const tmp<scalarField> tYOxidantpI =
//        patch().patchInternalField(YOxidant.internalField());

//    const scalarField& YOxidantpI = tYOxidantpI();
        
    // Return const access to internal cell values adjacent to patch of 
    // reductant mass fraction    
    const volScalarField& CReductant =
        db().lookupObject<volScalarField>(electrodeName_+"Reductant");
               

    const fvPatchField<scalar>& CReductantp =
        patch().patchField<volScalarField, scalar>(CReductant);
      
 
    //const tmp<scalarField> tYReductantpI =
    //    patch().patchInternalField(YReductant.internalField());

    //const scalarField& YReductantpI = tYReductantpI();
    

    const tmp<vectorField>& tn = patch().nf();
    vectorField n = tn();
//    const Field<scalar>& magS = patch().magSf();
           
    //vectorField jp(phiEp.size(), vector::zero);
    scalarField patchPhiE(phiEp.size(), 0.0);
    
    //***EDIT THIS INCLUDING NEWTON'S METHOD TO SOLVE BUTLER-VOLMER EQUATION****
    if (j.dimensions() == dimCurrent/dimArea)
    {
  
        const scalar tol = 1.e-6;
        const label iterMax = 100;        
        forAll(phiEp, faceI)
        {
            // Local electrode overpotential
            scalar overPot = electrodePotential_
                          -phiEp[faceI]-equilibriumPotential_;
            Info << "overPot_n = " << overPot << endl;                                                                  

            if (Tp[faceI] > 0.0)
            {
                scalar redFactor = pow(CReductantp[faceI]
                    /CRefReductant_,gammaReductant_);
                scalar oxFactor = pow(COxidantp[faceI]
                    /CRefOxidant_,gammaOxidant_);             

                if (redFactor > 0.0 || oxFactor > 0.0)
                {
                    scalar overPot_Old = overPot;                                  
                    scalar jps = jp[faceI] & n[faceI]*-1.0;
                    Info << "jps = " << jps << endl;                    
                    scalar K = F/R/Tp[faceI];                            
                    label iter = 0;
                    scalar eps = 1000.0;
                    scalar f_BV = 0.0;
                    scalar df_BV = 0.0;                                          
                    while ( eps > tol && iter < iterMax )
                    {           
                        //- Function of Butler-Volmer equation 
                        //  according to Newton's method            
//                        f_BV = jEx_*(redFactor*Foam::exp(alphaA_*K*overPot)
//                            -oxFactor*Foam::exp(-alphaC_*K*overPot))-jps; 
                            
                        f_BV = -jEx_*(Foam::exp(alphaA_*K*overPot)
                            -Foam::exp(-alphaC_*K*overPot))+jps;
                        //- Derivative of Bulter-Volmer equation
                        //  according to Newton's method                            
//                        df_BV = jEx_*K
//                            *(redFactor*alphaA_*Foam::exp(alphaA_*K*overPot)
//                            +oxFactor*alphaC_*Foam::exp(-alphaC_*K*overPot));
                            
                        df_BV = jEx_*K
                            *(alphaA_*Foam::exp(alphaA_*K*overPot)
                            +alphaC_*Foam::exp(-alphaC_*K*overPot));                            
                
                        overPot -= f_BV/df_BV;
                
                        eps = abs(overPot - overPot_Old);
                
                        overPot_Old = overPot;
                        iter++;
                    }
                    Info << "overPot_n+1 = " << overPot << endl;                                        
                    if (iter == iterMax-1)
                    {
                        Info << "Newton iterations did not converge within " << iter 
                        << " iterations \n for patch " << this->patch().name() 
                        << " of field " << this->dimensionedInternalField().name()
                        << " in file " << this->dimensionedInternalField().objectPath()
                        << " of face " << faceI << ".\n"
                        << " Residual error is " << eps;                
                    }                        
                }
            }
            else
            {
                FatalErrorIn("butlerVolmerPotentialFvPatchScalarField::updateCoeffs()")
                << "Temperature <= 0.0"
                << "\n    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
            }
            
            //overPot = 0.2;
            
              
            patchPhiE[faceI] = electrodePotential_
                -overPot-equilibriumPotential_;
        }                
        operator==(patchPhiE);
    }
    else
    {
        FatalErrorIn("butlerVolmerPotentialFvPatchScalarField::updateCoeffs()")
            << "dimensions of j are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::butlerVolmerPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "j", "j", jName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

/*
void Foam::butlerVolmerPotentialFvPatchScalarField::operator=
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
        butlerVolmerPotentialFvPatchScalarField
    );
}

// ************************************************************************* //

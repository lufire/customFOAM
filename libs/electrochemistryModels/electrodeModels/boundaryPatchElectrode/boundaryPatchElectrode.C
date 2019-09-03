/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "boundaryPatchElectrode.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

defineTypeNameAndDebug(boundaryPatchElectrode, 0);

addToRunTimeSelectionTable
(
    electrodeModel,
    boundaryPatchElectrode,
    electrolyte 
);

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

void boundaryPatchElectrode::initializeMembers()
{}

void boundaryPatchElectrode::correctElectrodeSurface()
{
    const labelUList& faceCells = 
        mesh_.boundary()[patchName_].faceCells();
    if (faceCells.size() == mesh_.boundary()[patchName_].magSf().size()) 
    {
        forAll(faceCells, faceI)
        {
            surfaceDensity_[faceCells[faceI]] = 
                mesh_.boundary()[patchName_].magSf()[faceI];
        }
    }
    else
    {
        scalar avgMagSf = 0.0;
        forAll(mesh_.boundary()[patchName_].magSf(), faceI)
        {
            avgMagSf += mesh_.boundary()[patchName_].magSf()[faceI];
        }
        avgMagSf /= mesh_.boundary()[patchName_].magSf().size(); 
        forAll(faceCells, faceI)
        {
            surfaceDensity_[faceCells[faceI]] = avgMagSf;
        }
    }
    surfaceDensity_ /= mesh_.V();
}

void boundaryPatchElectrode::correctElectrodeConductivity()
{}

void boundaryPatchElectrode::correctTransferCurrentDensity()
{
    if (CRef_[oxID_] <= 0.0 || CRef_[redID_] <= 0.0)
    {
        FatalErrorIn("boundaryPatchElectrode::calculateTransferCurrentDensity()")
            << "Reference concentrations must be greater than zero"
            << exit(FatalError);
    }
    else
    {
        // Physical constants
        const scalar& R = constant::physicoChemical::R.value();
        const scalar& F = constant::physicoChemical::F.value();

        //const scalarField& phiEls = resistanceModel_->phiEls();
        const scalarField& phiEls = 
            electrolyte_.phiE().boundaryField()[patchID_];

        //const scalarField& Coxs = resistanceModel_->Cs(oxID_);
        //const scalarField& Creds = resistanceModel_->Cs(redID_);
        const scalarField& Coxs = C_[oxID_].boundaryField()[patchID_];
        const scalarField& Creds = C_[redID_].boundaryField()[patchID_];

        Info << "Coxs: " << Coxs << endl;
        Info << "Creds: " << Creds << endl;

        const scalarField& T = 
            electrolyte_.thermo().T().boundaryField()[patchID_];

        vectorField& i = electrolyte_.i().boundaryField()[patchID_];

        const tmp<vectorField> tn = mesh_.boundary()[patchName_].nf();
        const vectorField& n = tn();



        scalar magSf = 0.0;
        scalar redFactor = 0.0;
        scalar oxFactor = 0.0;
        scalar i0 = 0.0;
        scalar iAvg = 0.0;
        scalar totalMagSf = 0.0;

        scalar eps = 1000;
        scalar tol = 1e-4;
        label iter = 0;
        label maxIter = 100;
        scalar f_BV = 0.0;
        scalar df_BV = 0.0;
        scalar eta = 0.0;
        scalar phiEs = avgPhiEs_;
        //const scalar eta_max = 0.5;
        const scalar minFactor = 0.1;
        //scalar phiEsOld = phiEsAvg_;

        while(eps > tol && iter < maxIter)
        {

            iAvg = 0.0;
            totalMagSf = 0.0;

            f_BV = 0.0;
            df_BV = 0.0;

            forAll(n, faceI)
            {
                
                magSf = mesh_.boundary()[patchName_].magSf()[faceI];
                totalMagSf += magSf; 
                
                scalar K = R*T[faceI]/(electronNumber_*F);
                //i[faceI] = currD_*-1.0*n[faceI];
                redFactor = 
                    pow(Creds[faceI]/CRef_[redID_],reactionOrder_[redID_]);
                oxFactor = 
                    pow(Coxs[faceI]/CRef_[oxID_],reactionOrder_[oxID_]);
                
                if(redFactor<minFactor)
                {
                    redFactor = minFactor;
                    //f_BV += -currD_*magSf;
                    //df_BV += 0.0;
                    //i[faceI] = n[faceI]*0.0;
                    //Info << "Concentrations too low: Coxs = " << Coxs[faceI] << ", Creds = " << Creds[faceI] << endl;
                }
                if(oxFactor<minFactor)
                {
                    oxFactor = minFactor;
                    //f_BV += -currD_*magSf;
                    //df_BV += 0.0;
                    //i[faceI] = n[faceI]*0.0;
                    //Info << "Concentrations too low: Coxs = " << Coxs[faceI] << ", Creds = " << Creds[faceI] << endl;
                }
                    eta_[faceI] = phiEs - phiEls[faceI] - eqPotential_
                        + K*log(oxFactor/redFactor);
                    Info << "faceI: " << faceI << ", eta: " << eta_[faceI] << ", phiEs: " << phiEs << ", phiEls: " << phiEls[faceI] << ", oxFactor: " << oxFactor << ", redFactor: " << redFactor << endl;
                    //Info << "faceI: " << faceI << ", eta: " << eta_[faceI] << endl;
                    //if(eta_[faceI] > 2.0)
                    //    eta_[faceI] = 2.0;
                    //if(eta_[faceI] < -2.0)
                    //    eta_[faceI] = -2.0;
                    //Info << "eta_[" << faceI << "] = " << eta_[faceI] << endl; 
                        
                    i0 = iEx_*pow(oxFactor,1.0-alphaA_)*pow(redFactor,alphaA_);

                    scalar iT = i0 
                       *(
                           exp(alphaA_/K*eta_[faceI])
                           -exp(-(1-alphaA_)/K*eta_[faceI])
                        );
                    f_BV += (iT - currD_)*magSf;
                    df_BV += (i0/K*(alphaA_*Foam::exp(alphaA_/K*eta_[faceI])
                        + (1.0-alphaA_)*Foam::exp(-(1.0-alphaA_)/K*eta_[faceI])))*magSf;

                    //d_eta = f_BV/df_BV;
                    //f_BV += f_BV*magSf;
                    //df_BV += df_BV*magSf;

                    //Info << "i0["<<faceI<<"] = " << i0 << endl;
                    i[faceI] = n[faceI]*-1.0*iT;
                    //scalar K = electronNumber_*F/(R*T[faceI]);
                        //eta_[faceI] = log(currD_/i0)/(alphaA_*K);
                    //i[faceI] = i0*-1.0*n[faceI]*Foam::exp(alphaA_*K*eta_[faceI]);p
                    iAvg += magSf*iT;
                    //Info << "i["<<faceI<<"] = " << i[faceI] << endl;
            }
            //if(totalMagSf != 0.0)
            //{

                //Info << "totalMagSf = " << totalMagSf << endl;
                f_BV /= totalMagSf;
                df_BV /= totalMagSf;
                iAvg /= totalMagSf;
            //}

            
            if(df_BV != 0.0)
            {
                phiEs -= f_BV/df_BV;
            }
            Info << "currD: " << currD_ << endl;
            Info << "iAvg: " << iAvg << endl;
            eps = abs(iAvg - currD_)/currD_;
            Info << "eps: " << eps << endl;
            iter++;
        }
        if (iter >= maxIter-1)
        {
            Info << "Newton iterations did not converge within " << iter 
            << " iterations \n for electrode of patch " 
            << mesh_.boundary()[patchName_].name()
            << " Residual error is " << eps << endl;                
        }                        
        //if(phiEs > 1.0)
        //    phiEs = 1.0;
        //if(phiEs < -1.0)
        //    phiEs = -1.0;
        avgPhiEs_ = phiEs; 
        forAll(phiEs_,faceI)
        {
            phiEs_ = phiEs;
        }
    }
}


void boundaryPatchElectrode::correctPotential()
{

    //iTrans_ *= 0.0;

    if (CRef_[oxID_] <= 0.0 || CRef_[redID_] <= 0.0)
    {
        FatalErrorIn("boundaryPatchElectrode::calculateTransferCurrentDensity()")
            << "Reference concentrations must be greater than zero"
            << exit(FatalError);
    }
    else
    {
        // Physical constants
        const scalar& R = constant::physicoChemical::R.value();
        const scalar& F = constant::physicoChemical::F.value();

        const scalarField& phiEls = resistanceModel_->phiEls();

        //const scalarField& Coxs = resistanceModel_->Cs(oxID_);
        //const scalarField& Creds = resistanceModel_->Cs(redID_);
        const scalarField& Coxs = C_[oxID_].boundaryField()[patchID_];
        const scalarField& Creds = C_[redID_].boundaryField()[patchID_];

        const scalarField& T = 
            electrolyte_.thermo().T().boundaryField()[patchID_];

        vectorField& i = electrolyte_.i().boundaryField()[patchID_];

        const tmp<vectorField> tn = mesh_.boundary()[patchName_].nf();
        const vectorField& n = tn();

        const scalar tol = 1.0e-6;
        const label iterMax = 100;
        const label& ne = electronNumber_;

        const scalar alphaC = 1.0-alphaA_;

        //scalarField eta = phiEs_;
        scalar etaGuess;
        //scalar i1;
        //scalar redFactor1;
        //scalar oxFactor1;
        //Info << "Cox = " << Coxs[5] << endl;
        //Info << "Cred = " << Creds[5] << endl;

        forAll(n, faceI)
        {

            //Info << "Creds["<<faceI<<"] = " << Creds[faceI] << endl;
            //Info << "Coxs["<<faceI<<"] = " << Coxs[faceI] << endl;
            scalar redFactor = 
                pow(Creds[faceI]/CRef_[redID_],reactionOrder_[redID_]);
            scalar oxFactor = 
                pow(Coxs[faceI]/CRef_[oxID_],reactionOrder_[oxID_]);

            if(redFactor <= 0.0) 
            {
                redFactor = SMALL;
            }
            if(oxFactor <= 0.0)
            {
                oxFactor = SMALL;
            }

            if(T[faceI] > 0.0)
            {
                //eta_[faceI] = phiEs_[faceI] - phiEls[faceI] - eqPotential_
                //    + R*T[faceI]/(ne*F)*log(redFactor/oxFactor);
                scalar i0 = 
                    iEx_*pow(oxFactor,alphaC)*pow(redFactor,alphaA_);
                scalar K = ne*F/(R*T[faceI]);

                //- First guess of overpotential through Tafel equation
                eta_[faceI] = log(currD_/i0)/(alphaA_*K);

                //if(faceI == 5)
                //{
                //    etaGuess = eta_[faceI];
                //    i1 = i0;
                //    redFactor1 = redFactor;
                //    oxFactor1 = oxFactor;
                //}


                ////scalar is = i[faceI] & n[faceI]*-1.0;
                //label iter = 0;
                //scalar eps = 1000.0;
                //scalar f_BV = 0.0;
                //scalar df_BV = 0.0;
                //scalar etaOld = eta_[faceI];

                //while(eps > tol && iter < iterMax)
                //{
                //    f_BV = i0*(Foam::exp(alphaA_*K*eta_[faceI])
                //        - Foam::exp(-alphaC*K*eta_[faceI]))
                //        - currD_;
                //    df_BV = i0*K*(alphaA_*Foam::exp(alphaA_*K*eta_[faceI])
                //        + alphaC*Foam::exp(-alphaC*K*eta_[faceI]));

                //    eta_[faceI] -= f_BV/df_BV;

                //    eps = abs(eta_[faceI] - etaOld);
                //    etaOld = eta_[faceI];
                //    iter++;
                //}
                //if (iter == iterMax-1)
                //{
                //    Info << "Newton iterations did not converge within " << iter 
                //    << " iterations \n for electrode of patch " 
                //    << mesh_.boundary()[patchName_].name()
                //    << " on face " << faceI << ".\n"
                //    << " Residual error is " << eps;                
                //}                        
            }
            else
            {
            FatalErrorIn("boundaryPatchElectrode::correctPotential()")
                << "Temperature must be > 0.0 "
                << "\n on patch " 
                << mesh_.boundary()[patchName_].name()
                << " on face " << faceI << ".\n"
                << exit(FatalError);
            }
            phiEs_[faceI] = eta_[faceI] + phiEls[faceI] + eqPotential_
                - R*T[faceI]/(ne*F)*log(oxFactor/redFactor);

        }

        // Calculate average electrode potential
        avgPhiEs_ = 0.0;
        forAll(phiEs_,faceI)
        {
            avgPhiEs_ += phiEs_[faceI];
        }
        if(phiEs_.size() > 0)
        {
            avgPhiEs_ /= phiEs_.size();
        }
    
        // Assign average electrode potential to field 
        forAll(phiEs_,faceI)
        {
            phiEs_[faceI] = avgPhiEs_;
        }

        correctTransferCurrentDensity();

        //Info << "i0 = " << i1 << endl;
        //Info << "oxFactor = " << oxFactor1 << endl;
        //Info << "redFactor = " << redFactor1 << endl;
        //Info << "etaGuess = " << etaGuess << endl;
        //Info << "eta_ = " << eta_[5] << endl;
        //Info << "phiEls_ = " << phiEls[5] << endl;
        //Info << "phiEs_ = " << phiEs_[5] << endl;
        //Info << "phiEls_ = " << phiEls[0] << endl;


        
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryPatchElectrode::boundaryPatchElectrode
(
    const dictionary& dict, 
    electrolyteModels::electrolyteModel& electrolyte 
)
:
    electrodeModel(dict, electrolyte),
    phiEs_(mesh_.boundary()[patchName_].magSf()),
    eta_(mesh_.boundary()[patchName_].magSf()),
    avgPhiEs_(readScalar(dict.lookup("potential"))),
    currD_(readScalar(dict.lookup("currentDensity"))),
    resistanceModel_(),
    fluxCorrector_(false)
{   
    phiEs_ /= mesh_.boundary()[patchName_].magSf();
    phiEs_ *= avgPhiEs_;
    eta_ *= 0.0;

    resistanceModel_.set
    (
        resistanceModel::New
        (
            dict.subDict("resistanceModel"),
            *this
        ).ptr()
    );
    correctElectrodeSurface();

}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

boundaryPatchElectrode::~boundaryPatchElectrode()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void boundaryPatchElectrode::correctElectricity()
{
    resistanceModel_->correctPotential();
    //correctElectrodeSurface();
    //correctElectrodeConductivity();
    if(currentControl_ == true)
    {
        correctPotential();
    }
    else
    {
        correctTransferCurrentDensity();
    }
    //correctSpeciesFlux();
}

void boundaryPatchElectrode::correctSpeciesFlux()
{
    const scalar& F = constant::physicoChemical::F.value();
    PtrList<volVectorField>& N = electrolyte_.N();
    //const volVectorField& i = electrolyte_.i();
    const vectorField& i =  electrolyte_.i().boundaryField()[patchID_];

    forAll(N,specieI)
    {
        vectorField& Ni = N[specieI].boundaryField()[patchID_];
        Ni = -i*stoichCoeff_[specieI]/(electronNumber_*F);
        //if(specieI == 0 && fluxCorrector_ == true)
        //{
        //    resistanceModel_->correctSpeciesFlux(Ni);
        //}
    }


}

void boundaryPatchElectrode::correctSpecies(PtrList<volScalarField>& C)
{
    const tmp<vectorField> tn = mesh_.boundary()[patchName_].nf();
    const vectorField& n = tn(); 

    const scalar Cstd = electrolyte_.Cstd().value();
    // Hard coded precipitation reaction
    forAll(C,specieI)
    {
        scalarField& Cp = C[specieI].boundaryField()[patchID_];
        forAll(Cp, faceI)
        {
            if(Cp[faceI] < Cstd*1e-3)
            {
                Cp[faceI] = Cstd*1e-3;
            }
        }
    }
    resistanceModel_->correctSpecies();

    //bool fluxCorrector = true;
    if(fluxCorrector_ == false)
    {
        scalar sCrit = 3.0;
        const scalarField& Cref = C[2].boundaryField()[patchID_];
        scalarField Csat = Cref;

        scalarField& Cp = C[0].boundaryField()[patchID_];
        //scalarField& Cs = resistanceModel_->Cs(oxID_);
        PtrList<scalarField>& Cs = resistanceModel_->Cs();
        forAll(Csat,faceI)
        {
            if(Cref[faceI]> 2.1*Cstd)
            {
                Csat[faceI] = -0.21*Cstd + 0.975e-1*Cref[faceI]
                    + 0.125e-2*(sqr(Cref[faceI])/Cstd);
            }
            else
            {
                Csat[faceI] = 20.0;
            }
            if(Csat[faceI]>1000)
            {
                Csat[faceI] = 1000;
            }

            if(Cp[faceI] > sCrit*Csat[faceI])
            {
                Cp[faceI] = sCrit*Csat[faceI];
            }
            //Info << "Cp["<<faceI<<"] = " << Cp[faceI] << endl;
            //Info << "Csat["<<faceI<<"] = " << Csat[faceI] << endl;
        }
        forAll(Cs,specieI)
        {
            forAll(Cs[specieI],faceI)
            {
                if(Cs[specieI][faceI] > sCrit*Csat[faceI])
                {
                    Cs[specieI][faceI] = sCrit*Csat[faceI];
                }
            }
        }
    }
}

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam

// ************************************************************************* //

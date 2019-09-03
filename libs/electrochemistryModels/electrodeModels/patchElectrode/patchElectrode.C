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

#include "patchElectrode.H"
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

defineTypeNameAndDebug(patchElectrode, 0);

addToRunTimeSelectionTable
(
    electrodeModel,
    patchElectrode,
    electrolyte 
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchElectrode::patchElectrode
(
    const dictionary& dict, 
    electrolyteModels::electrolyteModel& electrolyte
)
:
    electrodeModel(dict, electrolyte),
    phiEs_(readScalar(dict.lookup("phiEs")))
{
    //volScalarField surfaceDensity
    //(
    //    IOobject
    //    (
    //        "surfaceDensity",
    //        mesh_.time().timeName(),
    //        mesh_,
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //    ),
    //    mesh_,
    //    dimensionedScalar("surfaceDensity", dimArea, 0.0)
    //);
    
    labelUList faceCells = 
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


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

patchElectrode::~patchElectrode()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void patchElectrode::correctTransferCurrentDensity()
{
    const volScalarField& phiEl = electrolyte_.phiE();

    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();

    iTrans_ *= 0.0;

    if (CRef_[oxID_] <= 0.0 || CRef_[redID_] <= 0.0)
    {
        FatalErrorIn("patchElectrode::calculateTransferCurrentDensity()")
            << "Reference concentrations must be greater than zero"
            << exit(FatalError);
    }
    else
    {
        if (phiEl.dimensions() == dimMass*dimArea/pow3(dimTime)/dimCurrent)
        {
            labelUList faceCells = 
                mesh_.boundary()[patchName_].faceCells();
            label cellI;
            forAll(faceCells, faceI)
            {
                cellI = faceCells[faceI];
            /*
                Info << "jp = " << jp[faceI]
                << ", n = " << n[faceI]
                << ", T = " << Tp[faceI]            
                << ", CReductantp = " << CReductantp[faceI]                       
                << ", COxidantp = " << COxidantp[faceI]
                << ", CRefReductant = " << CRefReductant_
                << ", CRefOxidant = " << CRefOxidant_ << endl;
            */
                iTrans_[cellI] += surfaceDensity_[cellI]*iEx_
                   *(
                        pow(C_[redID_][cellI]/CRef_[redID_],
                            reactionOrder_[redID_])
                       *exp(alphaA_*F/R/T_[cellI]
                       *(phiEs_.value()-phiEl[cellI]-eqPotential_))
                       -pow(C_[oxID_][cellI]/CRef_[oxID_],
                           reactionOrder_[oxID_])
                       *exp(-alphaC_*F/R/T_[cellI]
                       *(phiEs_.value()-phiEl[cellI]-eqPotential_))
                    );
                //jp[faceI] =  n[faceI]*pow(YReductantp[faceI]/CRefReductant_,stoichCoeffReductant_);
                            
            }
            //Info << "iTrans = " << iTrans_[0] << endl; 
            //Info << "surfaceDensity = " << surfaceDensity_[0] << endl; 
            
        }
        else
        {
            FatalErrorIn(
                    "patchElectrode::calculateTransferCurrentDensity()")
                << "dimensions of field "
                << phiEl.name() 
                << "are not correct"
                << exit(FatalError);
        }
    }
    iTrans_.correctBoundaryConditions();
}

void patchElectrode::correctElectrodeSurface()
{}

void patchElectrode::correctElectrodeConductivity()
{}

void patchElectrode::correct()
{
    //correctElectrodeSurface();
    //correctElectrodeConductivity();
    correctTransferCurrentDensity();
}
//{
//    return true;
//}

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam

// ************************************************************************* //

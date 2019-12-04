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

#include "porousElectrode.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "fvm.H"
#include "myCfdemCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrodeModels
{
    defineTypeNameAndDebug(porousElectrode, 0);

    addToRunTimeSelectionTable
    (
        electrodeModel,
        porousElectrode,
        dictionary 
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electrodeModels::porousElectrode::porousElectrode
(
    const dictionary& dict, 
    electrolyteModels::electrolyteModel& electrolyte,
    const myCfdemCloud& particleCloud
)
:
    electrodeModel(dict, electrolyte),
    particleCloud_(particleCloud),
    percolationField_
    (
        IOobject
        (
            "percolation",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("percolation", dimless, 0.0)
    ),
    phiEs_
    (
        IOobject
        (
            "phiEs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    percolationInterval_(readScalar
    (
        electrodeDict_.lookup("percolationInterval")
    )),
    intervalCounter_(0),
    contactDistance_(readScalar
    (
        electrodeDict_.lookup("contactDistance")
    ))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::electrodeModels::porousElectrode::~porousElectrode()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
             
void Foam::electrodeModels::porousElectrode::correctTransferCurrentDensity
(
   const volScalarField& phiEl 
)
{
    
    // Physical constants
    const scalar& R = constant::physicoChemical::R.value();
    const scalar& F = constant::physicoChemical::F.value();

    iTrans_ *= 0.0;

    if (CRef_[oxID_] <= 0.0 || CRef_[redID_] <= 0.0)
    {
        FatalErrorIn("porousElectrode::calculateTransferCurrentDensity()")
            << "Reference concentrations must be greater than zero"
            << exit(FatalError);
    }
    else
    {
        if (phiEl.dimensions() == dimMass*dimArea/pow3(dimTime)/dimCurrent)
        {
            
            //Info << "CRefReductant = " << CRef_[redID_] << endl;
            //Info << "CRefOxidant = " << CRef_[oxID_] << endl;
            //Info << "T = " << T_[0] << endl;
            forAll(iTrans_, cellI)
            {
                if(surfaceDensity_[cellI] > SMALL)
                {
                
                    scalar potentialDiff = phiEs_[cellI] - phiEl[cellI];
                    //if (abs(potentialDiff) < abs(eqPotential_))
                    //{
                    iTrans_[cellI] = surfaceDensity_[cellI]*iEx_
                       *(
                            pow(C_[redID_][cellI]/CRef_[redID_],
                                reactionOrder_[redID_])
                           *Foam::exp(alphaA_*F/R/T_[cellI]
                           *(potentialDiff-eqPotential_))
                           -pow(C_[oxID_][cellI]/CRef_[oxID_],
                                reactionOrder_[oxID_])
                           *Foam::exp(-alphaC_*F/R/T_[cellI]
                           *(potentialDiff-eqPotential_))
                        );
                }
                else
                {
                    iTrans_[cellI] = 0.0;
                }
                //}
                //jp[faceI] =  n[faceI]*pow(YReductantp[faceI]/CRefReductant_,stoichCoeffReductant_);
                            
            }
            //iTrans_ = iTrans_*(0.0001) + iTrans_.oldTime()*(1-0.0001);
            //Info << "iTrans_[0] = " << iTrans_[0] << endl; 
            //Info << "surfaceDensity_[0] = " << surfaceDensity_[0] << endl; 
            //Info << "iTrans_[1] = " << iTrans_[1] << endl; 
            //Info << "surfaceDensity_[1] = " << surfaceDensity_[1] << endl; 
            //Info << "iTrans_[2] = " << iTrans_[2] << endl; 
            //Info << "surfaceDensity_[2] = " << surfaceDensity_[2] << endl; 
        }
        else
        {
            FatalErrorIn(
                    "porousElectrode::calculateTransferCurrentDensity()")
                << "dimensions of field "
                << phiEl.name() 
                << "are not correct"
                << exit(FatalError);
        }
    }
    iTrans_.correctBoundaryConditions();
}

void Foam::electrodeModels::porousElectrode::correctPercolation()
{
    if(intervalCounter_ == percolationInterval_-1)
    {
        label nParticles = particleCloud_.numberOfParticles();
        std::vector<bool> wallContact(nParticles, false);
        std::vector<label> particleContact(2);
        std::vector<std::vector<label> > particleContacts;
        scalar distance;

        // Populate wall contact list
        const fvPatch& patch = mesh_.boundary()[patchName_]; 
        //Info << "patchName: " << patchName_ << endl;
        //Info << "number of faces: " << patch.Cf().size() << endl;
        //Info << "patchName: " << "wall" << endl;
        //Info << "number of faces: " << mesh_.boundary()["wall"].Cf().size() << endl;
        //Info << "patchName: " << "inlet" << endl;
        //Info << "number of faces: " << mesh_.boundary()["inlet"].Cf().size() << endl;

        labelUList faceCells = patch.faceCells();
        label cellID;
        //Info << "faceCells: " << faceCells << endl;
        for(label index=0; index<nParticles; ++index)
        {
            cellID = particleCloud_.cellIDs()[index][0];
            if (cellID > -1)
            {
                //Info << "cellID: " << cellID << endl;
                forAll(faceCells, faceI)
                {
                    if(cellID == faceCells[faceI])
                    {
                        distance = ((-1.0*patch.nf()()[faceI])
                           &(particleCloud_.position(index) - patch.Cf()[faceI]))
                           /mag(patch.nf()()[faceI]);
                        if (distance <= particleCloud_.radius(index) + contactDistance_)
                        {
                            wallContact[index] = true;
                        }
                    }
                }
            }
        }
        // Populate inter-particle contact list
        for(label i=0; i<nParticles; ++i)
        {
            for(label j=0; j<nParticles; ++j)
            {
                if(i != j)
                {
                    distance = 
                       mag(particleCloud_.position(i) - particleCloud_.position(j));
                    if (distance <= (particleCloud_.radius(i)
                        +particleCloud_.radius(j) + contactDistance_))
                    {
                        particleContact[0] = i;
                        particleContact[1] = j;
                        particleContacts.push_back(particleContact);
                    }                
                }
            }
        }
        // Check indirect contacts to wall
        bool newContact = true;
        while (newContact)
        {
            newContact = false;
            for(label index=0; index<particleContacts.size(); ++index)
            {
                label i = particleContacts[index][0];
                label j = particleContacts[index][1];
                if( (wallContact[i] && !wallContact[j]) ||
                    (!wallContact[i] && wallContact[j]))
                {
                    wallContact[i] = true;
                    wallContact[j] = true;
                    newContact = true;
                }
            }
        }
        
        // Calculate percolation field
        percolationField_ *= 0.0;
        volScalarField particleSurfaceField(percolationField_);
        scalar particleSurface;
        for(label index=0; index<nParticles; ++index)
        {
            cellID = particleCloud_.cellIDs()[index][0];
            if(cellID > -1)
            {   
                particleSurface
                    = 4.0*constant::mathematical::pi
                    *pow(particleCloud_.radius(index),2.0);
                particleSurfaceField[cellID] 
                    += particleSurface;
                if (wallContact[index])
                {
                    percolationField_[cellID] 
                        += particleSurface;
                }
            }
        }
        forAll(percolationField_, cellI)
        {
            if(particleSurfaceField[cellI] > 0.0)
            {
                percolationField_[cellI] /= particleSurfaceField[cellI];
            }
        }
        intervalCounter_ = 0;
    }
    else
    {
        intervalCounter_++;
    }
}

void Foam::electrodeModels::porousElectrode::correctElectrodeSurface()
{
    //scalar totalParticleSurface(0.0);
    //scalar totalCellVolume(0.0);
    //scalarField undividedVoidfraction(mesh_.V());
    scalar particleSurface(0.0);
    surfaceDensity_ *= 0.0;
    label cellID = -1;
    for(label index=0; index<particleCloud_.numberOfParticles(); ++index)
    {
        cellID = particleCloud_.cellIDs()[index][0];
        if(cellID > -1)
        {
            particleSurface = 4.0*constant::mathematical::pi
                *pow(particleCloud_.radius(index),2.0);
            //totalParticleSurface += particleSurface; 
            surfaceDensity_[cellID] += particleSurface;
            //totalCellVolume += mesh_.V()[cellID];
            //undividedVoidfraction[cellID] -= 4.0/3.0
            //    *constant::mathematical::pi*pow(particleRadius,3.0);
        }
    }
    surfaceDensity_ /= mesh_.V();
    surfaceDensity_ *= percolationField_.internalField();
    //surfaceDensity_ = 100.0;

    //Info << "Total Specific Electrode Surface = " << totalParticleSurface/totalCellVolume << endl; 
    //undividedVoidfraction /= mesh_.V();
    //surfaceDensity_ *= 
    //    particleCloud_.voidFractionM().voidFractionNext().internalField()
    //    /undividedVoidfraction;
}

void Foam::electrodeModels::porousElectrode::correctElectrodeConductivity()
{
    sigmaEff_ = sigma_*pow(1.0-particleCloud_.voidFractionM().voidFractionNext(),1.5);
    //sigmaEff_ = sigma_*pow(1.0-0.8,1.5);
}

void Foam::electrodeModels::porousElectrode::correct(const volScalarField& phiEl)
{
    correctPercolation();
    correctElectrodeSurface();
    correctElectrodeConductivity();
    correctTransferCurrentDensity(phiEl);

    // Solve for electric potential of electrode
    fvScalarMatrix phiEsEqn
    (
        fvm::laplacian(sigmaEff_, phiEs_, "laplacian(sigma,phiE)") == iTrans_
    );
    //phiEs *= 0;
    phiEsEqn.solve();
    phiEs_.correctBoundaryConditions();
    // Calculate solid current density
    is_ = -sigmaEff_*fvc::grad(phiEs_);
    is_.correctBoundaryConditions();
}
//bool Foam::porousElectrode::read()
//{
//    return true;
//}


// ************************************************************************* //

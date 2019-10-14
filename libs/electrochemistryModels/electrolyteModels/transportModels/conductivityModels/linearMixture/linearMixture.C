/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "linearMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "electrolyteModel.H"
//#include "diffusivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(linearMixture, 0);
addToRunTimeSelectionTable
(
    conductivityModel,
    linearMixture, 
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearMixture::linearMixture
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    conductivityModel(dict, electrolyte),
    balanceIndex_(0),
    z_(electrolyte.z()),
    C_(electrolyte.thermo().composition().C())
{
    label nIons = 0;
    forAll(z_,i)
    {
        if(z_[i] != 0)
        {
            nIons++;
        }
    }
    //- Last charged species is assumed to balance previous charged species 
    balanceIndex_ = nIons-1;

    if(balanceIndex_ <= 0)
    {
        FatalErrorIn
        (
            "conductivityModel::linearMixture::linearMixture("
            "const dictionary&, const electrolyteModel&)"
        )   << "Number of charged species (ions) must be at least 2"
            << exit(FatalError);
    }

    lambda_.setSize(nIons-1);
    word balanceIon = species_[balanceIndex_];
    forAll(lambda_,i)
    {
       
        word couplename = balanceIon + word("-") + species_[i];
        if
        (
            dict.subDict("equivalentConductivity")
                .lookupEntryPtr(couplename,1,1) == NULL
        )
        {
            couplename = species_[i] + word("-") + balanceIon;
        }
        
        if(species_[i] != balanceIon)
        {
            lambda_.set
            (
                i, 
                new dimensionedScalar
                (
                    dict.subDict("equivalentConductivity")
                        .lookup(couplename)
                )
            );
        }

    }

    X_.setSize(species_.size());
    forAll(X_,i)
    {
        X_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "X_" + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("X", dimless, 0.0)
            )
        );
    }

    update();
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void linearMixture::update()
{
    volScalarField rX(C_[0]*abs(z_[0]));
    volScalarField oneC 
    (                
        IOobject
        (
            "oneC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("oneC", dimMoles/dimVol, 1.0)
    );

    for(label i = 1; i<X_.size(); ++i)
    {
        if(i != balanceIndex_ && z_[i] != 0)
        {
            rX += C_[i]*abs(z_[i]);
        }
    }
    forAll(X_,i)
    {
        if(i == balanceIndex_)
        {
            X_[i] = 1.0;
        }
        else if(z_[i] != 0)
        {
            X_[i] = (C_[i]+oneC*SMALL)*abs(z_[i])/(rX+oneC*SMALL); 
        }
        else
        {
            X_[i] = 0.0;
        }
    }

    volScalarField lambda(X_[0]*lambda_[0]);
    for(label i = 1; i<lambda_.size(); ++i)
    {
        lambda += X_[i]*lambda_[i];
    }
    kappa_ = lambda*(C_[balanceIndex_]+oneC*SMALL);

    // Set internal field next to boundary equal to boundary values
    // (leads to better results for weighted potential gradient calculation)
    forAll(kappa_.boundaryField(), patchI)
    {
        fvPatchScalarField& kappap = kappa_.boundaryField()[patchI];
        const fvPatch& patch = kappa_.boundaryField()[patchI].patch();
        const labelUList& faceCells = patch.faceCells();
        forAll(kappap, faceI)
        {
            kappa_[faceCells[faceI]] = kappap[faceI];
        }
    }
}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

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

#include "speciesConductivity.H"
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

defineTypeNameAndDebug(speciesConductivity, 0);
addToRunTimeSelectionTable
(
    transferenceNumberModel,
    speciesConductivity, 
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

speciesConductivity::speciesConductivity
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    transferenceNumberModel(dict, electrolyte),
    balanceIndex_(0),
    z_(electrolyte.z()),
    C_(electrolyte.thermo().composition().C())
{
    speciesKappa_.setSize(species_.size());
    forAll(speciesKappa_,i)
    {
        if(z_[i] != 0)
        {
            speciesKappa_.set
            (
                i, 
                new dimensionedScalar
                (
                    dict.subDict("speciesConductivity")
                        .lookup(species_[i])
                )
            );
        }
        else
        {
            speciesKappa_.set
            (
                i, 
                new dimensionedScalar
                (
                    "zero",
                    sqr(dimCurrent)*pow3(dimTime)/(dimMoles*dimMass),
                    0.0
                )
            );
        }

    }
    update();
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void speciesConductivity::update()
{
    volScalarField one 
    (                
        IOobject
        (
            "one",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", pow3(dimTime)/dimVol/dimMass*sqr(dimCurrent), 1.0)
    );
    volScalarField rt(C_[0]*abs(z_[0])*speciesKappa_[0]);
    for(label i = 1; i<C_.size(); ++i)
    {
        rt += C_[i]*abs(z_[i])*speciesKappa_[i];
    }

    forAll(t_,i)
    {
        t_[i] = C_[i]*abs(z_[i])*speciesKappa_[i]/(rt+one*SMALL);
    }
}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

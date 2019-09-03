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

#include "constantDiff.H"
#include "electrolyteModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels 
{

defineTypeNameAndDebug(constantDiff, 0);
addToRunTimeSelectionTable
(
    diffusivityModel,
    constantDiff, 
    constant 
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantDiff::constantDiff
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    diffusivityModel(dict, thermo)
{
    forAll(species_,speciesI)
    {
        D_[speciesI] = dict.subDict("ionicConducitivity")
            .lookup(species_[speciesI]);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

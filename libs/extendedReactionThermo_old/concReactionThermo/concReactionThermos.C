/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "makeReactionThermo.H"

#include "concReactionThermo.H"
#include "heRhoThermo.H"

#include "specie.H"
#include "hConstThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"

#include "constTransport.H"

#include "concMultiComponentMixture.H"
#include "concReactingMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Multi-component thermo for internal energy
//makeReactionMixtureThermo
//(
//    rhoThermo,
//    concReactionThermo,
//    heRhoThermo,
//    concMultiComponentMixture,
//    constEThermoPhysics
//);

//// Multi-component reaction thermo
//makeReactionMixtureThermo
//(
//    rhoThermo,
//    concReactionThermo,
//    heRhoThermo,
//    concReactingMixture,
//    constEThermoPhysics
//);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    concReactionThermo,
    heRhoThermo,
    concReactingMixture,
    constEThermoPhysics
);

// Multi-component thermo for sensible enthalpy
//makeReactionMixtureThermo
//(
//    rhoThermo,
//    concReactionThermo,
//    heRhoThermo,
//    concMultiComponentMixture,
//    constHThermoPhysics
//);

//// Multi-component reaction thermo
//makeReactionMixtureThermo
//(
//    rhoThermo,
//    concReactionThermo,
//    heRhoThermo,
//    concReactingMixture,
//    constHThermoPhysics
//);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    concReactionThermo,
    heRhoThermo,
    concReactingMixture,
    constHThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

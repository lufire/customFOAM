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

#ifndef makeChemistrySolverTypes_H
#define makeChemistrySolverTypes_H

#include "chemistrySolver.H"

#include "extendedChemistryModel.H"

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeConcChemistrySolverType(SS, Comp, Thermo)                             \
                                                                              \
    typedef SS<extendedChemistryModel<Comp, Thermo> > SS##Comp##Thermo;               \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##Thermo,                                                     \
        (#SS"<" + word(Comp::typeName_())                                     \
      + "," + Thermo::typeName() + ">").c_str(),                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Comp,                                                                 \
        SS##Comp##Thermo,                                                     \
        fvMesh                                                                \
    );


#define makeConcChemistrySolverTypes(CompChemModel,Thermo)                        \
                                                                              \
    makeConcChemistrySolverType                                                   \
    (                                                                         \
        noChemistrySolver,                                                    \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \
                                                                              \
    makeConcChemistrySolverType                                                   \
    (                                                                         \
        EulerImplicit,                                                        \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \
                                                                              \
    makeConcChemistrySolverType                                                   \
    (                                                                         \
        ode,                                                                  \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::liquidSchmidtNumber

Description
    Multispecies mass transport model based on Schmidt number.
  
SourceFiles
    liquidSchmidtNumber.C

\*---------------------------------------------------------------------------*/

#ifndef liquidSchmidtNumber_H
#define liquidSchmidtNumber_H

#include "liquidFick.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class liquidSchmidtNumber Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class liquidSchmidtNumber
:
    public liquidFick<ThermoType>
{  

protected:
  
    // Protected data

        //- Schmidt number
        scalar Sc_;

    // Protected member functions
    
        //- Update the coefficients of the species equations
        void updateCoefficients();
    

public:

    //- Runtime type information
    TypeName("liquidSchmidtNumber");

    
    // Constructors

        //- Construct from components
        liquidSchmidtNumber
        (
            // A.Alexiou 2014
            //hsCombustionThermo& thermo,
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        );


    //- Destructor
    ~liquidSchmidtNumber()
    {}
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "liquidSchmidtNumber.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

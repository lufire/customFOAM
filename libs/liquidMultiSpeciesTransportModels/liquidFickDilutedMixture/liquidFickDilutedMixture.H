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
    Foam::liquidFickDilutedMixture

Description
    Fick multispecies mass transport model for very diluted mixture.

SourceFiles
    liquidFickDilutedMixture.C

\*---------------------------------------------------------------------------*/

#ifndef liquidFickDilutedMixture_H
#define liquidFickDilutedMixture_H

#include "liquidFick.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class Bosanquet Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class liquidFickDilutedMixture
:
    public liquidFick<ThermoType>
{

protected:

    // Protected member functions
    
        //- Update the coefficients of the species equations
        //void updateCoefficients();


public:

    //- Runtime type information
    TypeName("liquidFickDilutedMixture");

    
    // Constructors

        //- Construct from components
        liquidFickDilutedMixture
        (
            // A.Alexiou 2014
            //hsCombustionThermo& thermo,
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        );


    //- Destructor
    ~liquidFickDilutedMixture()
    {}
    
    // Member Functions
    
        //- Update the coefficients of the species equations
        void updateCoefficients();    
        
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "liquidFickDilutedMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

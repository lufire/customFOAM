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
    Foam::liquidBosanquet

Description
    Fick multispecies mass transport model with Knudsen diffusivity correction.
  
SourceFiles
    liquidBosanquet.C

\*---------------------------------------------------------------------------*/

#ifndef liquidBosanquet_H
#define liquidBosanquet_H

#include "liquidFick.H"
#include "KnudsenDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class liquidBosanquet Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class liquidBosanquet
:
    public liquidFick<ThermoType>
{

protected:
  
    // Protected data
    
        //- Diffusivity models 
        autoPtr<KnudsenDiffusivityModel> DKModel_;


    // Protected member functions
    
        //- Update the coefficients of the species equations
        void updateCoefficients();
  
public:

    //- Runtime type information
    TypeName("liquidBosanquet");

    
    // Constructors

        //- Construct from components
        liquidBosanquet
        (
            // A.Alexiou 2014
            // hsCombustionThermo& thermo,
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        );


    //- Destructor
    ~liquidBosanquet()
    {}
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "liquidBosanquet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


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
    Foam::concentratedElectrolyte

Description
    concentratedElectrolyte multispecies mass transport model.
  
SourceFiles
    concentratedElectrolyte.C

\*---------------------------------------------------------------------------*/

#ifndef concentratedElectrolyte_H
#define concentratedElectrolyte_H

#include "electrolyteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels 
{

/*---------------------------------------------------------------------------*\
                             Class concentratedElectrolyte Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class concentratedElectrolyte
:
    public electrolyteModel
{  

protected:
  
    // Protected data
    
        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;

    // Protected member functions
    
        //- Update the coefficients of the species equations
        virtual void updateCoefficients();
    

public:

    //- Runtime type information
    TypeName("concentratedElectrolyte");

    
    // Constructors

        //- Construct from components
        concentratedElectrolyte
        (
            // A.Alexiou 2014
            // hsCombustionThermo& thermo,
            const dictionary& dict,
            const concReactionThermo& thermo
        );

        //- Destructor
        ~concentratedElectrolyte()
        {}

    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depend on thermoType template
        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }
    
        //- Specific enthalpy for a specie given by index and temperature.
        //  Defined here becouse depend on thermoType template
        scalar hs(label i, label T) const
        {
            return speciesThermo_[i].hs(i,T);
        }
        
        //- Update electrolyte properties 
        virtual void correct();   

        //- Update the coefficients of the species equations
        //virtual void updateCoefficients();            

        //- Solve the mass transport equantions without reactions
        //virtual scalar correct
        //(
        //    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        //);
                       
        //- Read transportProperties dictionary
        //virtual bool read();        
    
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "concentratedElectrolyte.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

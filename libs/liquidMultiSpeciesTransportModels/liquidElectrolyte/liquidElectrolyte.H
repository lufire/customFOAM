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
    Foam::liquidElectrolyte

Description
    liquidElectrolyte multispecies mass transport model.
  
SourceFiles
    liquidElectrolyte.C

\*---------------------------------------------------------------------------*/

#ifndef liquidElectrolyte_H
#define liquidElectrolyte_H

#include "fvc.H"
#include "fvm.H"
#include "liquidMultiSpeciesTransportModel.H"
#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class liquidElectrolyte Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class liquidElectrolyte
:
    public liquidMultiSpeciesTransportModel
{  

protected:
  
    // Protected data
    
        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;

        //- Ionic conductivity model 
        autoPtr<conductivityModel> kappaModel_;

        //- liquidMultiSpecies diffusion coefficients
        autoPtr<volScalarField> kappa_;

    // Protected member functions
    
        //- Update the coefficients of the species equations
//        virtual void updateCoefficients();
    

public:

    //- Runtime type information
    TypeName("liquidElectrolyte");

    
    // Constructors

        //- Construct from components
        liquidElectrolyte
        (
            // A.Alexiou 2014
            // hsCombustionThermo& thermo,
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        );


      //- Destructor
      ~liquidElectrolyte()
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
        
        //- Update the coefficients of the species equations
        virtual void updateCoefficients();            

        //- Solve the mass transport equantions without reactions
        virtual scalar correct
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );
                       
        //- Read transportProperties dictionary
        virtual bool read();        
    
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "liquidElectrolyte.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

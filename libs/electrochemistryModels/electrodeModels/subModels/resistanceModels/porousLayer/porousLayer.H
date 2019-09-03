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
    Foam::resistanceModel::porousLayer

Description
    Constant diffusivity model.

SourceFiles
    porousLayer.C

\*---------------------------------------------------------------------------*/

#ifndef porousLayer_H
#define porousLayer_H

#include "resistanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels 
{
/*---------------------------------------------------------------------------*\
                           Class porousLayer Declaration
\*---------------------------------------------------------------------------*/

class porousLayer
:
    public resistanceModel
{
    // Private data

        //- Porosity of layer
        scalar eps_;

        //- Thickness of layer
        scalar delta_;

        //- Thickness of layer
        scalarField deltaField_;

        //- Bruggemann coefficient
        scalar brugg_;


    // Private Member Functions
        
        //- 
    
public:

    //- Runtime type information
    TypeName("porousLayer");


    // Constructors

        //- construct from components
        porousLayer
        (
            const dictionary& dict,
            const electrodeModel& electrode 
        );


    // Destructor

        ~porousLayer()
        {}


    // Member Functions

        //- Correct species concentrations 
        virtual void correctSpecies();

        //- Correct electric potential 
        virtual void correctPotential();

        //- Correct the species flux 
        void correctSpeciesFlux(vectorField& N);

        //- Update surface values (concentrations and electric potential) 
        //inline void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrodeModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// *********************************************************************** //

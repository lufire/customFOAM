/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more dnuils.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Potanin

Description
    Thixotropic flow model according to Potanin (2010).

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    Potanin.C

\*---------------------------------------------------------------------------*/

#ifndef Potanin_H
#define Potanin_H

//#include "thixotropyModel.H"
#include "PotaninStructure.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Potanin Declaration
\*---------------------------------------------------------------------------*/

class Potanin
:
    public rheologyModel
{

    // Private data

        // Model constants

        //- Parameter for viscosity function 
        const dimensionedScalar tau0_;

        //- Parameter for viscosity function 
        const dimensionedScalar sr0_;

        //- Dimensionless exponent viscosity function
        const dimensionedScalar k_;

        //- Reference to Potanin structure model
        PotaninStructure& PotaninStructure_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Potanin(const Potanin&);

        //- Disallow default bitwise assignment
        void operator=(const Potanin&);

public:

    //- Runtime type information
    TypeName("Potanin");


    // Constructors

        //- Construct from components
        Potanin
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Potanin()
        {}


    // Member Functions
        
        //- Correct the rheology stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
#endif

// ************************************************************************* //

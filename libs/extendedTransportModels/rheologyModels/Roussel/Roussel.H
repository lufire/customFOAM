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
    Roussel

Description
    Thixotropic flow model according to Roussel (2010).

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    Roussel.C

\*---------------------------------------------------------------------------*/

#ifndef Roussel_H
#define Roussel_H

//#include "thixotropyModel.H"
//#include "RousselStructure.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Roussel Declaration
\*---------------------------------------------------------------------------*/

class Roussel
:
    public rheologyModel
{

    // Private data

        // Model constants

        //- Yield stress equivalent to the value in the equilibrium function
        const dimensionedScalar tau0_;

        //- Regularization parameter according to Papanastasiou
        const dimensionedScalar m_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Roussel(const Roussel&);

        //- Disallow default bitwise assignment
        void operator=(const Roussel&);

public:

    //- Runtime type information
    TypeName("Roussel");


    // Constructors

        //- Construct from components
        Roussel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Roussel()
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

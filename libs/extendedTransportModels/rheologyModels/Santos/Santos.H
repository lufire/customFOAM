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
    Santos

Description
    Thixotropic flow model according to Santos et al. (2017).

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    Santos.C

\*---------------------------------------------------------------------------*/

#ifndef Santos_H
#define Santos_H

//#include "thixotropyModel.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Santos Declaration
\*---------------------------------------------------------------------------*/

class Santos
:
    public rheologyModel
{

    // Private data

        // Model constants

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Santos(const Santos&);

        //- Disallow default bitwise assignment
        void operator=(const Santos&);

public:

    //- Runtime type information
    TypeName("Santos");


    // Constructors

        //- Construct from components
        Santos
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Santos()
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

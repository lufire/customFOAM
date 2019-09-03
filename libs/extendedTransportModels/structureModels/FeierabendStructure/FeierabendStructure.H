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
    FeierabendStructure

Description
    Thixotropic structure model after Andrei FeierabendStructure 
    (A. FeierabendStructure, 2010)

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    FeierabendStructure.C

\*---------------------------------------------------------------------------*/

#ifndef FeierabendStructure_H
#define FeierabendStructure_H

#include "structureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FeierabendStructure Declaration
\*---------------------------------------------------------------------------*/

class FeierabendStructure
:
    public structureModel
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        // Model constants

        //- Fully structured state
        const dimensionedScalar lambda0_;

        //- Characteristic time for structure build-up 
        const dimensionedScalar theta_;

        //- Parameter for structure build-up 
        const dimensionedScalar a_;

        //- Parameter for structure build-up 
        const dimensionedScalar b_;

        //- Parameter for structure break-down 
        const dimensionedScalar alpha_;

        //- Parameter for structure break-down 
        const dimensionedScalar c_;

        //- Parameter for structure break-down 
        const dimensionedScalar d_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        FeierabendStructure(const FeierabendStructure&);

        //- Disallow default bitwise assignment
        void operator=(const FeierabendStructure&);


public:

    //- Runtime type information
    TypeName("FeierabendStructure");

    // Constructors

        //- Construct from components
        FeierabendStructure
        (
            const word& name,
            const dictionary& dict,
            const rheologyModel& rheology,
            const volVectorField& U
        );


    // Destructor

        virtual ~FeierabendStructure()
        {}


    // Member Functions
    //
        //- Correct the structure field 
        virtual void correct(const surfaceScalarField& phi);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

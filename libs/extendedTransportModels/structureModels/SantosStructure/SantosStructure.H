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
    SantosStructure

Description
    Thixotropic structure model after Andrei SantosStructure 
    (A. SantosStructure, 2010)

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    SantosStructure.C

\*---------------------------------------------------------------------------*/

#ifndef SantosStructure_H
#define SantosStructure_H

#include "structureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SantosStructure Declaration
\*---------------------------------------------------------------------------*/

class SantosStructure
:
    public structureModel
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        // Model constants

        //- Factor for break-down/build-up rate
        const dimensionedScalar k1_;

        //- Exponent for break-down/build-up rate
        const dimensionedScalar n1_;

        //- Factor for equilibrium structure parameter
        const dimensionedScalar k2_;

        //- Exponent for equilibrium structure parameter
        const dimensionedScalar n2_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SantosStructure(const SantosStructure&);

        //- Disallow default bitwise assignment
        void operator=(const SantosStructure&);


public:

    //- Runtime type information
    TypeName("SantosStructure");

    // Constructors

        //- Construct from components
        SantosStructure
        (
            const word& name,
            const dictionary& dict,
            const rheologyModel& rheology,
            const volVectorField& U
        );


    // Destructor

        virtual ~SantosStructure()
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

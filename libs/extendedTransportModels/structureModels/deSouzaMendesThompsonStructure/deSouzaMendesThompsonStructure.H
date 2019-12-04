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
    deSouzaMendesThompsonStructure

Description
    Thixotropic structure model after de Souza Mendes and Thompson
    (Paulo R. de Souza Mendes & Roney L. Thompson, 2013)

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    deSouzaMendesThompsonStructure.C

\*---------------------------------------------------------------------------*/

#ifndef deSouzaMendesThompsonStructure_H
#define deSouzaMendesThompsonStructure_H

#include "structureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class deSouzaMendesThompsonStructure Declaration
\*---------------------------------------------------------------------------*/

class deSouzaMendesThompsonStructure
:
    public structureModel
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        // Model constants

        //- Relaxation time 
        const dimensionedScalar teq_;

        //- Model exponent 
        const dimensionedScalar a_;

        //- Model exponent 
        const dimensionedScalar b_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        deSouzaMendesThompsonStructure(const deSouzaMendesThompsonStructure&);

        //- Disallow default bitwise assignment
        void operator=(const deSouzaMendesThompsonStructure&);


public:

    //- Runtime type information
    TypeName("deSouzaMendesThompsonStructure");

    // Constructors

        //- Construct from components
        deSouzaMendesThompsonStructure
        (
            const word& name,
            const dictionary& dict,
            const rheologyModel& rheology,
            const volVectorField& U
        );


    // Destructor

        virtual ~deSouzaMendesThompsonStructure()
        {}


    // Member Functions

        //- Correct the structure field 
        virtual void correct(const surfaceScalarField& phi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

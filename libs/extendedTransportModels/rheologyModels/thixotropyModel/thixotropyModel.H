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
    thixotropyModel

Description
    Abstract base class for the thixotropic flow model.

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    thixotropyModel.C

\*---------------------------------------------------------------------------*/

#ifndef thixotropyModel_H
#define thixotropyModel_H

#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class thixotropyModel Declaration
\*---------------------------------------------------------------------------*/

class thixotropyModel
:
    public rheologyModel 
{
protected:
    // Protected data

        //- Kinetic structure model
        autoPtr<structureModel> structureModelPtr_;

        // Model constants

        // Fields

private:
    // Private Member Functions

        //- Disallow default bitwise copy construct
        thixotropyModel(const thixotropyModel&);

        //- Disallow default bitwise assignment
        void operator=(const thixotropyModel&);

public:

    //- Runtime type information
    TypeName("thixotropyModel");


    // Constructors

        //- Construct from components
        thixotropyModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~thixotropyModel()
        {}

        
        //- Correct the rheology stress
        virtual void correct() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
#endif

// ************************************************************************* //

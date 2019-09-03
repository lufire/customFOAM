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
    viscoElasticModel

Description
    Abstract base class for viscoelastic flow models.

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    viscoElasticModel.C

\*---------------------------------------------------------------------------*/

#ifndef viscoElasticModel_H
#define viscoElasticModel_H

#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoElasticModel Declaration
\*---------------------------------------------------------------------------*/

class viscoElasticModel
:
    public rheologyModel 
{
protected:
    // Protected data

        // Model constants

        // Fields

        //- Transported viscoElastic stress
        volSymmTensorField tau_;

private:
    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoElasticModel(const viscoElasticModel&);

        //- Disallow default bitwise assignment
        void operator=(const viscoElasticModel&);

public:

    //- Runtime type information
    TypeName("viscoElasticModel");


    // Constructors

        //- Construct from components
        viscoElasticModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~viscoElasticModel()
        {}


    // Member Functions
        
        //- Return the rheology stress tensor
        inline tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Correct the rheology stress
        virtual void correct() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//
#endif

// ************************************************************************* //

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
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    viscoElasticModel

Description
    Base class for viscoElastic flow models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    viscoElasticModel.C
    newViscoelasticModel.C

\*---------------------------------------------------------------------------*/

#ifndef viscoElasticModel_H
#define viscoElasticModel_H

#include "viscoElasticLaw.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoElasticModel Declaration
\*---------------------------------------------------------------------------*/

class viscoElasticModel
:
    public IOdictionary
{
private:

    // Private data

        //- Viscoelastic law
        autoPtr<viscoElasticLaw> lawPtr_;


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

        //- Return the viscoElastic stress tensor
        virtual tmp<volSymmTensorField> tau() const;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

        //- Correct the viscoElastic stress
        virtual void correct();

        //- Read viscoElasticProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

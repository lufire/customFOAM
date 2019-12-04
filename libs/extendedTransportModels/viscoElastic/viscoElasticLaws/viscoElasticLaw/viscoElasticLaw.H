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
    viscoElasticLaw

Description
    Abstract base class for viscoElastic flow models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    viscoElasticLaw.C
    newViscoelasticLaw.C

\*---------------------------------------------------------------------------*/

#ifndef viscoElasticLaw_H
#define viscoElasticLaw_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoElasticLaw Declaration
\*---------------------------------------------------------------------------*/

class viscoElasticLaw
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to face flux field
        const surfaceScalarField& phi_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoElasticLaw(const viscoElasticLaw&);

        //- Disallow default bitwise assignment
        void operator=(const viscoElasticLaw&);


public:

    //- Runtime type information
    TypeName("viscoElasticLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscoElasticLaw,
            dictionary,
            (
                const word& name,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const dictionary& dict
            ),
            (name, U, phi, dict)
        );


    // Constructors

        //- Construct from components
        viscoElasticLaw
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Selectors

        //- Return a reference to the selected viscoElastic law
        static autoPtr<viscoElasticLaw> New
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~viscoElasticLaw()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the face flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return the viscoElastic stress tensor
        virtual tmp<volSymmTensorField> tau() const = 0;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const = 0;

        //- Correct the viscoElastic stress
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

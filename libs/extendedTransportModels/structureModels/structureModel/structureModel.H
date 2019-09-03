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
    structureModel

Description
    Base class for thixotropic structure models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    structureModel.C
    newStructureModel.C

\*---------------------------------------------------------------------------*/

#ifndef structureModel_H
#define structureModel_H

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

//Forward declaration of classes
class rheologyModel;

/*---------------------------------------------------------------------------*\
                        Class structureModel Declaration
\*---------------------------------------------------------------------------*/

class structureModel
{
protected:

    // Private data

        //- Name
        const word name_;

        //- Model subdictionary 
        const dictionary dict_;

        //- Reference to parent class
        const rheologyModel& rheology_;

        //- Structure parameter field
        volScalarField lambda_;

        //- Explicit structure parameter source
        volScalarField expLambdaSource_;

        //- Implicit structure parameter source
        fvScalarMatrix impLambdaSource_;

private:
    // Private Member Functions

        //- Disallow default bitwise copy construct
        structureModel(const structureModel&);

        //- Disallow default bitwise assignment
        void operator=(const structureModel&);


public:

    //- Runtime type information
    TypeName("structureModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            structureModel,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const rheologyModel& rheology,
                const volVectorField& U
            ),
            (name, dict, rheology, U)
        );


    // Constructors

        //- Construct from components
        structureModel
        (
            const word& name,
            const rheologyModel& rheology,
            const volVectorField& U
        );


    // Selectors

        //- Return a reference to the selected structure model 
        static autoPtr<structureModel> New
        (
            const word& name,
            const rheologyModel& rheology,
            const volVectorField& U
        );


    // Destructor

        virtual ~structureModel()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the structure field
        const volScalarField& lambda() const
        {
            return lambda_;
        }

        //- Correct the viscoElastic stress
        virtual void correct(const surfaceScalarField& phi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

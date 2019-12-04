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
    PotaninStructure

Description
    Thixotropic structure model after Andrei PotaninStructure 
    (A. PotaninStructure, 2010)

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    PotaninStructure.C

\*---------------------------------------------------------------------------*/

#ifndef PotaninStructure_H
#define PotaninStructure_H

#include "structureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PotaninStructure Declaration
\*---------------------------------------------------------------------------*/

class PotaninStructure
:
    public structureModel
{
    // Private data

        //- Dictionary
        const dictionary& dict_;

        // Model constants

        //- Recovery time constant 
        const dimensionedScalar tR_;

        //- Breakup time constant 
        const dimensionedScalar tB_;

        //- Reference shear rate 
        const dimensionedScalar srRef_;

        //- Model exponent for recovery source 
        const dimensionedScalar nR_;

        //- Model exponent for breakup sink 
        const dimensionedScalar nB_;

        //- Model exponent for characteristic recovery time
        const dimensionedScalar mR_;

        //- Model exponent for characteristic breakup time
        const dimensionedScalar mB_;

        // Fields

        //- Shear dependent factor field f
        volScalarField f_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PotaninStructure(const PotaninStructure&);

        //- Disallow default bitwise assignment
        void operator=(const PotaninStructure&);


public:

    //- Runtime type information
    TypeName("PotaninStructure");

    // Constructors

        //- Construct from components
        PotaninStructure
        (
            const word& name,
            const dictionary& dict,
            const rheologyModel& rheology,
            const volVectorField& U
        );


    // Destructor

        virtual ~PotaninStructure()
        {}


    // Member Functions
        
        //- Access to shear rate function f
        inline volScalarField& f()
        {
            return f_;
        }

        //- Correct the structure field 
        virtual void correct(const surfaceScalarField& phi);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

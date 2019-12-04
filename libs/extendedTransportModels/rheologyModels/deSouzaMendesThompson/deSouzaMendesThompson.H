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
    deSouzaMendesThompson

Description
    Base class for the rheology flow model.

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    deSouzaMendesThompson.C
    newViscoelasticModel.C

\*---------------------------------------------------------------------------*/

#ifndef deSouzaMendesThompson_H
#define deSouzaMendesThompson_H

//#include "thixoElastoViscoModel.H"
//#include "viscoElasticModel.H"
//#include "thixotropyModel.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class deSouzaMendesThompson Declaration
\*---------------------------------------------------------------------------*/

class deSouzaMendesThompson
:
    public rheologyModel 
    //public thixotropyModel
{

    // Private data

        // Model constants

        //- Elastic modulus at fully structured state
        const dimensionedScalar G0_;

        //- Dimensionless positive constant for elastic modulus calculation
        const dimensionedScalar m_;


        // Fields

        //- Shear stress field
        volSymmTensorField tau_;

        //- Shear rate tensor
        volSymmTensorField twoD_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        deSouzaMendesThompson(const deSouzaMendesThompson&);

        //- Disallow default bitwise assignment
        void operator=(const deSouzaMendesThompson&);

public:

    //- Runtime type information
    TypeName("deSouzaMendesThompson");


    // Constructors

        //- Construct from components
        deSouzaMendesThompson
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~deSouzaMendesThompson()
        {}


    // Member Functions
      
        //- Calculate divergence of shear stress field    
        tmp<fvVectorMatrix> divTau(const volVectorField& U) const;

        //- Return the rheology stress tensor
        inline tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Correct the rheology stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
#endif

// ************************************************************************* //

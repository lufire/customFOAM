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
    WhiteMetznerCross

Description
    White-Metzner Power Law non linear viscoElastic fluid model
    (Larson, R. G., Constitutive equations for polymer melts
    and solutions, 1988).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    WhiteMetznerCross.C

\*---------------------------------------------------------------------------*/

#ifndef WhiteMetznerCross_H
#define WhiteMetznerCross_H

#include "viscoElasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WhiteMetznerCross Declaration
\*---------------------------------------------------------------------------*/

class WhiteMetznerCross
:
    public viscoElasticLaw
{
    // Private data

        //- Transported viscoElastic stress
        volSymmTensorField tau_;


        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- Relaxation time
            dimensionedScalar lambda_;

            //- Fitted parameter for etaP
            dimensionedScalar K_;

            //- Fitted parameter for lambda
            dimensionedScalar L_;

            //- Fitted parameter for etaP
            dimensionedScalar m_;

            //- Fitted parameter for lambda
            dimensionedScalar n_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WhiteMetznerCross(const WhiteMetznerCross&);

        //- Disallow default bitwise assignment
        void operator=(const WhiteMetznerCross&);


public:

    //- Runtime type information
    TypeName("WhiteMetznerCross");

    // Constructors

        //- Construct from components
        WhiteMetznerCross
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~WhiteMetznerCross()
        {}


    // Member Functions

        //- Return the viscoElastic stress tensor
        virtual tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

        //- Correct the viscoElastic stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

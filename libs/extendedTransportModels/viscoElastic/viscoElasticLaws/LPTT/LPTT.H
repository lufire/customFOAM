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
    LPTT

Description
    Phan-Thin-Tanner (linear form) viscoElastic fluid model.
    (Phan-Thien, N.; Tanner, R. I., 1977).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    LPTT.C

\*---------------------------------------------------------------------------*/

#ifndef LPTT_H
#define LPTT_H

#include "viscoElasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LPTT Declaration
\*---------------------------------------------------------------------------*/

class LPTT
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

            //- Extensibility parameter
            dimensionedScalar epsilon_;

            //- Relaxation time
            dimensionedScalar lambda_;

            //- Chains slip in relation to the continue environment
            dimensionedScalar zeta_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LPTT(const LPTT&);

        //- Disallow default bitwise assignment
        void operator=(const LPTT&);


public:

    //- Runtime type information
    TypeName("PTT-Linear");

    // Constructors

        //- Construct from components
        LPTT
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~LPTT()
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

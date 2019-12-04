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
    rheologyModel

Description
    Abstract base class for all rheology flow models.

Author
    Lukas Feierabend.  All rights reserved

SourceFiles
    rheologyModel.C
    newRheologyodel.C

\*---------------------------------------------------------------------------*/

#ifndef rheologyModel_H
#define rheologyModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
//#include "singlePhaseTransportModel.H"
#include "viscosityModel.H"
#include "structureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class viscosityModel;

/*---------------------------------------------------------------------------*\
                        Class rheologyModel Declaration
\*---------------------------------------------------------------------------*/

class rheologyModel
:
    public IOdictionary
{

protected:

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to face flux field
        const surfaceScalarField& phi_;

        //- Viscosity model 
        autoPtr<viscosityModel> viscosityModelPtr_;

        // Model constants

        //- Density
        const dimensionedScalar rho_;

        //- Zero shear rate equilibrium viscosity
        const dimensionedScalar nu0_;

        //- Infinite shear rate equilibrium viscosity
        const dimensionedScalar nuInf_;

        //- Kinetic structure model
        autoPtr<structureModel> structureModelPtr_;

        // Fields

        //- Equilibrium viscosity
        const volScalarField& nuEq_;

        //- Shear rate 
        volScalarField sr_;

        //- Apparent viscosity
        volScalarField nu_;

    // Protected Member Functions
    //
private:

    // Private Member Functions
        
        //- Calculate strain rate field
        //tmp<volScalarField> strainRate();

        //- Disallow default bitwise copy construct
        rheologyModel(const rheologyModel&);

        //- Disallow default bitwise assignment
        void operator=(const rheologyModel&);

public:

    //- Runtime type information
    TypeName("rheologyModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rheologyModel,
            dictionary,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );

    // Constructors

        //- Construct from components
        rheologyModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    // Selectors

        //- Return a reference to the selected structure model 
        static autoPtr<rheologyModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~rheologyModel()
        {}

    // Member Functions

        //- Calculate divergence of shear stress field    
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U) const;

        //- Return apparent viscosity field
        inline volScalarField& nu();

        //- Return apparent viscosity field
        inline const volScalarField& nu() const;

        //- Return equilibrium viscosity field
        inline const volScalarField& nuEq() const;

        //- Return strain rate field
        inline const volScalarField& sr() const;

        //- Return zero shear-rate viscosity 
        inline const dimensionedScalar& nu0() const;

        //- Return infinite shear-rate viscosity 
        inline const dimensionedScalar& nuInf() const;

        //- Correct the rheology stress
        virtual void correct();

        //- Read rheologyProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rheologyModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::Casson

Description
     Herschel-Bulkley non-Newtonian viscosity model.

SourceFiles
    Casson.C

\*---------------------------------------------------------------------------*/

#ifndef Casson_H
#define Casson_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Casson Declaration
\*---------------------------------------------------------------------------*/

class Casson
:
    public viscosityModel
{
    // Private data

        dictionary CassonCoeffs_;

        dimensionedScalar nu0_;
        dimensionedScalar nuInf_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
	    dimensionedScalar srMin_;

        volScalarField nu_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("Casson");


    // Constructors

        //- Construct from components
        Casson
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~Casson()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

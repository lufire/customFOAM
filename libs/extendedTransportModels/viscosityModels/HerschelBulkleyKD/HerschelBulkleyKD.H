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
    Foam::viscosityModels::HerschelBulkleyKD

Description
     Herschel-Bulkley non-Newtonian viscosity model.

SourceFiles
    HerschelBulkleyKD.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkleyKD_H
#define HerschelBulkleyKD_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HerschelBulkleyKD Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkleyKD
:
    public viscosityModel
{
    // Private data

        dictionary HerschelBulkleyKDCoeffs_;

        dimensionedScalar k_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
        dimensionedScalar nu0_;
	    dimensionedScalar m_;
	    scalar maxSolidFraction_;
	    scalar EinsteinCoeff_;
        word voidfractionName_;
        const volScalarField& voidfraction_;
        IOobject voidfractionHeader_;
        volScalarField nu_;
        //mutable tmp<volScalarField> kinVisc_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu();

public:

    //- Runtime type information
    TypeName("HerschelBulkleyKD");


    // Constructors

        //- Construct from components
        HerschelBulkleyKD
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~HerschelBulkleyKD()
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

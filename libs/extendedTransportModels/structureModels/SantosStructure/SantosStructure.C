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

\*---------------------------------------------------------------------------*/

#include "SantosStructure.H"
#include "addToRunTimeSelectionTable.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SantosStructure, 0);
    addToRunTimeSelectionTable(structureModel, SantosStructure, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SantosStructure::SantosStructure
(
    const word& name,
    const dictionary& dict,
    const rheologyModel& rheology,
    const volVectorField& U
)
:
    structureModel(name, rheology, U),
    dict_(dict),
    k1_(dict.lookup("k1")),
    n1_(dict.lookup("n1")),
    k2_(dict.lookup("k2")),
    n2_(dict.lookup("n2"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SantosStructure::correct(const Foam::surfaceScalarField& phi)
{
    lambda_.max(SMALL);
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    const volScalarField& sr = rheology_.sr();

    const volScalarField c = k1_*pow(sr*tone, n1_)*rtone; 
    const volScalarField lambdaEq = k2_*pow(sr*tone, n2_);

    expLambdaSource_ = c*(2.0*lambda_*lambdaEq - pow(lambdaEq, 2.0));
    impLambdaSource_ = fvm::Sp(c*lambda_, lambda_); 

    structureModel::correct(phi);

    lambda_.max(1.0);
}

// ************************************************************************* //

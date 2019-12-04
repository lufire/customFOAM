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

#include "CoussotStructure.H"
#include "addToRunTimeSelectionTable.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CoussotStructure, 0);
    addToRunTimeSelectionTable(structureModel, CoussotStructure, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoussotStructure::CoussotStructure
(
    const word& name,
    const dictionary& dict,
    const rheologyModel& rheology,
    const volVectorField& U
)
:
    structureModel(name, rheology, U),
    dict_(dict),
    theta_(dict.lookup("theta")),
    alpha_(dict.lookup("alpha"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CoussotStructure::correct(const Foam::surfaceScalarField& phi)
{
    lambda_.max(SMALL);

    const volScalarField& sr = rheology_.sr();

    expLambdaSource_ = 1.0/theta_; 
    impLambdaSource_ = fvm::Sp(alpha_*sr, lambda_);

    structureModel::correct(phi);

    lambda_ = max(lambda_, 0.0);
    //lambda_ = min(lambda_, 500.0);

    Info<< "lambda min/max = "
        << gMin(lambda_) << "/"
        << gMax(lambda_) << endl;
}

// ************************************************************************* //

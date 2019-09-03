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

#include "generalizedNewtonian.H"
#include "addToRunTimeSelectionTable.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalizedNewtonian, 0);
addToRunTimeSelectionTable(rheologyModel, generalizedNewtonian, dictionary);


// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

generalizedNewtonian::generalizedNewtonian
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    rheologyModel(U, phi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void generalizedNewtonian::correct()
{
    // Update rheology base model 
    rheologyModel::correct();
    nu_ = nuEq_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

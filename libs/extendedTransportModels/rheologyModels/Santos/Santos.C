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

#include "Santos.H"
#include "addToRunTimeSelectionTable.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Santos, 0);
addToRunTimeSelectionTable(rheologyModel, Santos, dictionary);


// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Santos::Santos
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    rheologyModel(U, phi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Santos::correct()
{
    // Update rheology base model 
    rheologyModel::correct();

    //Update structure parameter 
    structureModelPtr_->correct(phi_);

    // Retrieve structure field
    const volScalarField& lambda = structureModelPtr_->lambda();

    //Update apparent viscosity
    nu_ = lambda*nuEq_;
    nu_ = max(nu_, nuInf_);
    nu_ = min(nu_, nu0_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

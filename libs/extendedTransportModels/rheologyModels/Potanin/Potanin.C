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

#include "Potanin.H"
#include "addToRunTimeSelectionTable.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Potanin, 0);
addToRunTimeSelectionTable(rheologyModel, Potanin, dictionary);


// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Potanin::Potanin
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    rheologyModel(U, phi),
    //thixotropyModel(U, phi),
    tau0_(this->subDict((this->type()+"Coeffs")).lookup("tau0")),
    sr0_(this->subDict((this->type()+"Coeffs")).lookup("sr0")),
    k_(this->subDict((this->type()+"Coeffs")).lookup("k")),
    PotaninStructure_(dynamic_cast<PotaninStructure&>(structureModelPtr_()))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Potanin::correct()
{
    // Update rheology base model 
    rheologyModel::correct();

    //Update structure parameter 
    volScalarField& f = PotaninStructure_.f();
    f = nuInf_*(1.0 + pow(sr0_/sr_, k_));
    structureModelPtr_->correct(phi_);

    // Retrieve structure field
    const volScalarField& lambda = structureModelPtr_->lambda();

    //Update apparent viscosity
    nu_ = lambda*f;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

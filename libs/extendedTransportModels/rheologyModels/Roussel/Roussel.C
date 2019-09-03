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

#include "Roussel.H"
#include "addToRunTimeSelectionTable.H"
//#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Roussel, 0);
addToRunTimeSelectionTable(rheologyModel, Roussel, dictionary);


// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Roussel::Roussel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    rheologyModel(U, phi),
    //thixotropyModel(U, phi),
    tau0_(this->subDict((this->type()+"Coeffs")).lookup("tau0")),
    m_(this->subDict((this->type()+"Coeffs")).lookup("m"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Roussel::correct()
{
    // Update rheology base model 
    rheologyModel::correct();

    //Update structure parameter 
    structureModelPtr_->correct(phi_);

    // Retrieve structure field
    const volScalarField& lambda = structureModelPtr_->lambda();

    //Update apparent viscosity
    nu_ = nuEq_;
    dimensionedScalar tone("tone", dimTime, 1.0);
    nu_ += (lambda*tau0_*(1.0-exp(-m_*tone*sr_)))
        /(max(sr_, dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)));
    nu_ = max(nu_, nuInf_);
    nu_ = min(nu_, nu0_);

    Info<< "nu min/max = "
        << gMin(nu_) << "/"
        << gMax(nu_) << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

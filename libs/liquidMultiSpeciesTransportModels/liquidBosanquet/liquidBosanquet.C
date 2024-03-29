/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "liquidBosanquet.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::liquidBosanquet<ThermoType>::updateCoefficients()
{
    liquidFick<ThermoType>::updateCoefficients();
    
    forAll(this->D_, i)
    {
      this->D_.set
      (
          i,
          1/(1/this->D_[i] + 1/this->thermo_.rho()/DKModel_().D(i))
      );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::liquidBosanquet<ThermoType>::liquidBosanquet
(
    // A.Alexiou 2014
    // hsCombustionThermo& thermo,
    concReactionThermo& thermo,
    const incompressible::turbulenceModel& turbulence
)
:
    liquidFick<ThermoType>(thermo, turbulence)
{
    // Construct the Knudsen diffusivity model
    DKModel_.set(new KnudsenDiffusivityModel(thermo.T(), this->species()));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels 
{

inline porosityModelList&
electrolyteModel::pZones()
{
    return pZones_;
}

inline const concReactionThermo&
electrolyteModel::thermo() const
{
    return thermo_;
}

inline const speciesTable&
electrolyteModel::species() const
{
    return thermo_.composition().species();
}

inline const PtrList<volVectorField>& 
electrolyteModel::N() const
{
    return N_;
}

inline PtrList<volVectorField>& 
electrolyteModel::N()
{
    return N_;
}

inline const volVectorField&
electrolyteModel::N(const label i) const
{
    return N_[i];
}

inline volVectorField&
electrolyteModel::N(const label i)
{
    return N_[i];
}

inline const label&
electrolyteModel::z(const label i) const
{
    return z_[i];
}

inline const PtrList<label>& 
electrolyteModel::z() const
{
    return z_;
}

inline const word& 
electrolyteModel::balanceSpecies() const
{
    return balanceSpecies_;
}

inline const word& 
electrolyteModel::inertSpecies() const
{
    return inertSpecies_;
}

inline const volScalarField&
electrolyteModel::D(label i) const
{
    return diffModel().D(i);
}

inline const PtrList<volScalarField>& 
electrolyteModel::D() const
{
    return diffModel().D();
}

inline const volScalarField&
electrolyteModel::kappa() const
{
    return kappaModel().kappa();
}

inline const dimensionedScalar&
electrolyteModel::Cstd() const
{
    return Cstd_;
}

inline const dimensionedScalar&
electrolyteModel::epsilon() const
{
    return epsilon_;
}

inline const diffusivityModel& 
electrolyteModel::diffModel() const
{
    return diffModel_();
}

inline const conductivityModel& 
electrolyteModel::kappaModel() const
{
    return kappaModel_();
}

inline const volScalarField&
electrolyteModel::phiE() const 
{
    return phiE_;
}

inline volScalarField&
electrolyteModel::phiE() 
{
    return phiE_;
}

inline const volVectorField&
electrolyteModel::i() const 
{
    return i_;
}

inline volVectorField&
electrolyteModel::i() 
{
    return i_;
}

inline const volVectorField&
electrolyteModel::iD() const 
{
    return iD_;
}

inline volVectorField&
electrolyteModel::iD() 
{
    return iD_;
}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam
// ************************************************************************* //

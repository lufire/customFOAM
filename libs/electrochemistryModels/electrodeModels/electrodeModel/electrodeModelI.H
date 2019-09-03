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
namespace electrodeModels
{

inline const electrolyteModels::electrolyteModel&
electrodeModel::electrolyte() const
{
    return electrolyte_;
}

inline volScalarField&
electrodeModel::iTrans()
{
    return iTrans_;
}

inline volVectorField&
electrodeModel::is()
{
    return is_;
}

inline const PtrList<label>&
electrodeModel::stoichCoeff() const
{
    return stoichCoeff_;
}

inline label
electrodeModel::stoichCoeff(const label i) const
{
    return stoichCoeff_[i];
}

inline label
electrodeModel::electronNumber() const
{
    return electronNumber_;
}

inline label
electrodeModel::oxidantIndex() const
{
    return oxID_;
}

inline label
electrodeModel::reductantIndex() const
{
    return redID_;
}

inline label
electrodeModel::patchIndex() const
{
    return patchID_;
}

inline const word&
electrodeModel::patchName() const
{
    return patchName_;
}

inline const volScalarField&
electrodeModel::sigmaEff() const
{
    return sigmaEff_;
}

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam
// ************************************************************************* //

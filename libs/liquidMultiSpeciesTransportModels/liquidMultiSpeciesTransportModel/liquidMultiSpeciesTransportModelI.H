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

inline const Foam::volScalarField&
Foam::liquidMultiSpeciesTransportModel::Dij(label i, label j) const
{
    return DijModel_().D(i,j);
}

inline const Foam::volScalarField&
Foam::liquidMultiSpeciesTransportModel::Di(label i) const
{
    return D_[i];
}

inline const Foam::PtrList<Foam::volScalarField>& 
Foam::liquidMultiSpeciesTransportModel::D() const
{
    return D_;
}

inline Foam::porosityModelList&
Foam::liquidMultiSpeciesTransportModel::pZones()
{
    return pZones_;
}


const inline Foam::speciesTable&
Foam::liquidMultiSpeciesTransportModel::species() const
{
    return thermo_.composition().species();
}


inline const Foam::surfaceScalarField&
Foam::liquidMultiSpeciesTransportModel::n(const label i) const
{
    return n_[i];
}

inline const Foam::label&
Foam::liquidMultiSpeciesTransportModel::z(const label i) const
{
    return z_[i];
}

inline const Foam::PtrList<Foam::label>& 
Foam::liquidMultiSpeciesTransportModel::z() const
{
    return z_;
}

// ************************************************************************* //

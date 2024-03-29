/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

inline const Foam::speciesTable&
Foam::concBasicMultiComponentMixture::species() const
{
    return species_;
}


inline bool Foam::concBasicMultiComponentMixture::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}


inline bool Foam::concBasicMultiComponentMixture::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>& Foam::concBasicMultiComponentMixture::active() const
{
    return active_;
}


inline void Foam::concBasicMultiComponentMixture::setActive(label speciei)
{
    active_[speciei] = true;
}


inline void Foam::concBasicMultiComponentMixture::setInactive(label speciei)
{
    active_[speciei] = false;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::concBasicMultiComponentMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::concBasicMultiComponentMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::concBasicMultiComponentMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::concBasicMultiComponentMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::concBasicMultiComponentMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::concBasicMultiComponentMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}

inline Foam::PtrList<Foam::volScalarField>&
Foam::concBasicMultiComponentMixture::C()
{
    return C_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::concBasicMultiComponentMixture::C() const
{
    return C_;
}


inline Foam::volScalarField& Foam::concBasicMultiComponentMixture::C(const label i)
{
    return C_[i];
}


inline const Foam::volScalarField& Foam::concBasicMultiComponentMixture::C
(
    const label i
) const
{
    return C_[i];
}


inline Foam::volScalarField& Foam::concBasicMultiComponentMixture::C
(
    const word& specieName
)
{
    return C_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::concBasicMultiComponentMixture::C
(
    const word& specieName
) const
{
    return C_[species_[specieName]];
}

// ************************************************************************* //

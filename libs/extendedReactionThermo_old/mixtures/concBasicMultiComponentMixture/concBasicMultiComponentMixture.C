/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "concBasicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::concBasicMultiComponentMixture::concBasicMultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    species_(specieNames),
    Y_(species_.size()),
    C_(species_.size()),
    initialRho_
    (
        IOobject
        ( 
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    )
{
    forAll(species_, i)
    {
        IOobject header
        (
            "C_"+species_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );
        if (header.typeHeaderOk<volScalarField>(true))
        {
            C_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "C_"+species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField Cdefault
            (
                IOobject
                (
                    "Cdefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            C_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "C_"+species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Cdefault
                )
            );
        }
        Y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    species_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("Ydefault", dimless, 0.0)
            )
        );
    }
    // Do not enforce constraint of sum of mass fractions to equal 1 here
    // - not applicable to all models
}


// ************************************************************************* //

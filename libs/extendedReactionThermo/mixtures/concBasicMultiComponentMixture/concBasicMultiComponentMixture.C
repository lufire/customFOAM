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

#include "concBasicMultiComponentMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(concBasicMultiComponentMixture, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::concBasicMultiComponentMixture::concBasicMultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMixture(thermoDict, mesh, phaseName),
    species_(specieNames),
    active_(species_.size(), true),
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
    tmp<volScalarField> tCdefault;

    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName("C_"+species_[i], phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            C_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("C_"+species_[i], phaseName),
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
            // Read Ydefault if not already read
            if (!tCdefault.valid())
            {
                word CdefaultName(IOobject::groupName("Cdefault", phaseName));

                IOobject timeIO
                (
                    CdefaultName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject constantIO
                (
                    CdefaultName,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject time0IO
                (
                    CdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.typeHeaderOk<volScalarField>(true))
                {
                    tCdefault = new volScalarField(timeIO, mesh);
                }
                else if (constantIO.typeHeaderOk<volScalarField>(true))
                {
                    tCdefault = new volScalarField(constantIO, mesh);
                }
                else
                {
                    tCdefault = new volScalarField(time0IO, mesh);
                }
            }

            C_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("C_"+species_[i], phaseName),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tCdefault()
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
                    IOobject::groupName(species_[i], phaseName),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
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

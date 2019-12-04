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
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "structureModel.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::structureModel> Foam::structureModel::New
(
    const word& name,
    const rheologyModel& rheology,
    const volVectorField& U
)
{
    if(rheology.found("structureModel"))
    {
        word typeName = rheology.lookup("structureModel");
        const dictionary& dict = rheology.subDict(typeName+"Coeffs");

        Info<< "Selecting thixotropic structure model " << typeName << endl;

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(typeName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "structureModel::New(const word& name, const volVectorField&, "
                "const surfaceScalarField&)"
            )   << "Unknown structureModel type " << typeName
                << endl << endl
                << "Valid structureModel types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<structureModel>(cstrIter()(name, dict, rheology, U));
    }
    else
    {
        autoPtr<structureModel> structureModelPtr(new structureModel(name, rheology, U));
        return structureModelPtr;
    }
}


// ************************************************************************* //

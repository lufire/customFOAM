/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::concChemistryModel

Description
    Chemistry model for density-based thermodynamics

SourceFiles
    concChemistryModelI.H
    concChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef concChemistryModel_H
#define concChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "concReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class concChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class concChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        concChemistryModel(const concChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const concChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<concReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("conc");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        concChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        concChemistryModel(const fvMesh& mesh);


    //- Selector
    static autoPtr<concChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~concChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline concReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const concReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "concChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

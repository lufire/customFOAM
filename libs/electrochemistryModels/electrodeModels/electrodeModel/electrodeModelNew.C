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

#include "electrodeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

autoPtr<electrodeModel>
electrodeModel::New
(
    const dictionary& dict, 
    electrolyteModels::electrolyteModel& electrolyte 
)
{

    const word modelType(dict.lookup("typeName"));

    Info<< "Selecting electrode model " << modelType << endl;

    electrolyteConstructorTable::iterator cstrIter =
        electrolyteConstructorTablePtr_->find(modelType);

    if (cstrIter == electrolyteConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "electrodeModel::New(const dictionary&, "
            "const electrolyteModel&)"
        )   << "Unknown electrodeModel type "
            << modelType << endl << endl
            << "Valid  electrodeModels are : " << endl
            << electrolyteConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<electrodeModel>
        (cstrIter()(dict, electrolyte));
}

//autoPtr<electrodeModel>
//electrodeModel::New
//(
//    const dictionary& electrodeDict, 
//    const concReactionThermo& thermo
//)
//{
//
//    const word modelType(electrodeDict.lookup("typeName"));
//
//    Info<< "Selecting electrode model " << modelType << endl;
//
//    dictionaryConstructorTable::iterator cstrIter =
//        dictionaryConstructorTablePtr_->find(modelType);
//
//    if (cstrIter == dictionaryConstructorTablePtr_->end())
//    {
//        FatalErrorIn
//        (
//            "electrodeModel::New(const dictionary&, "
//            "const concReactionThermo&)"
//        )   << "Unknown electrodeModel type "
//            << modelType << endl << endl
//            << "Valid  electrodeModels are : " << endl
//            << dictionaryConstructorTablePtr_->toc()
//            << exit(FatalError);
//    }
//
//    return autoPtr<electrodeModel>
//        (cstrIter()(electrodeDict, thermo));
//}

autoPtr<electrodeModel>
electrodeModel::New
(
    const dictionary& dict,
    electrolyteModels::electrolyteModel& electrolyte,
    const myCfdemCloud& particleCloud
)
{
    const word modelType(dict.lookup("typeName"));

    Info<< "Selecting electrode model " << modelType << endl;

    myCfdemCloudConstructorTable::iterator cstrIter =
        myCfdemCloudConstructorTablePtr_->find(modelType);

    if (cstrIter == myCfdemCloudConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "electrodeModel::New(const dictionary&, "
            "const electrolyteModel&, const myCfdemCloud&)"
        )   << "Unknown electrodeModel type "
            << modelType << endl << endl
            << "Valid  electrodeModels are : " << endl
            << myCfdemCloudConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<electrodeModel>
        (cstrIter()(dict, electrolyte, particleCloud));
}

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam
// ************************************************************************* //

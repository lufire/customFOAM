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

#include "diffusivityModel.H"
//#include "conductivityModel.H"
#include "electrolyteModel.H"
//#include "volFields.H"

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(diffusivityModel, 0);
defineRunTimeSelectionTable(diffusivityModel, dictionary);
defineRunTimeSelectionTable(diffusivityModel, constant);
//defineRunTimeSelectionTable(diffusivityModel, Electrolyte);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
diffusivityModel::diffusivityModel
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    dict_(dict),
    mesh_(electrolyte.thermo().T().mesh()),
    thermo_(electrolyte.thermo()),
    species_(thermo_.composition().species()),
    D_(species_.size())
{
    forAll(species_, specieI)
    {
        D_.set
        (
            specieI, new volScalarField
            (
                IOobject
                (
                    "D_" + species_[specieI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "D",
                    dimArea/dimTime,
                    0.0
                )
            )
        );
    }
}

diffusivityModel::diffusivityModel
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    dict_(dict),
    mesh_(thermo.T().mesh()),
    thermo_(thermo),
    species_(thermo.composition().species()),
    D_(species_.size())
{
    forAll(species_, specieI)
    {
        D_.set
        (
            specieI, new volScalarField
            (
                IOobject
                (
                    "D_" + species_[specieI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "D",
                    dimArea/dimTime,
                    0.0
                )
            )
        );
    }
}

//diffusivityModel::diffusivityModel
//(
//    const dictionary& dict,
//    const concReactionThermo& thermo,
//    const dimensionedScalar& Cstd 
//)
//:
//    dict_(dict),
//    mesh_(thermo.T().mesh()),
//    thermo_(thermo),
//    species_(thermo.composition().species())
//{
//    forAll(species_, specieI)
//    {
//        D_.set
//        (
//            specieI, new volScalarField
//            (
//                IOobject
//                (
//                    "D_" + species_[specieI],
//                    mesh_.time().timeName(),
//                    mesh_,
//                    IOobject::NO_READ,
//                    IOobject::AUTO_WRITE
//                ),
//                mesh_,
//                dimensionedScalar
//                (
//                    "D",
//                    sqr(dimArea)/dimTime,
//                    0.0
//                )
//            )
//        );
//    }
//}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<diffusivityModel> diffusivityModel::New
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
{
    word modelType(dict.lookup("typeName"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "diffusivityModel::New(const dictionary&, "
            "const electrolyteModel&)"
        )   << "Unknown diffusivityModel type "
            << modelType << endl << endl
            << "Valid diffusivityModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<diffusivityModel>
        (cstrIter()(dict, electrolyte));
}

autoPtr<diffusivityModel> diffusivityModel::New
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
{
    word modelType(dict.lookup("typeName"));

    constantConstructorTable::iterator cstrIter =
        constantConstructorTablePtr_->find(modelType);

    if (cstrIter == constantConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "diffusivityModel::New(const dictionary&, "
            "const concReactionThermo&)"
        )   << "Unknown diffusivityModel type "
            << modelType << endl << endl
            << "Valid diffusivityModels are : " << endl
            << constantConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<diffusivityModel>
        (cstrIter()(dict, thermo));
}

//autoPtr<diffusivityModel> diffusivityModel::New
//(
//    const dictionary& dict,
//    const concReactionThermo& thermo,
//    const dimensionedScalar& Cstd 
//)
//{
//    word modelType(dict.lookup("typeName"));
//
//    dictionaryConstructorTable::iterator cstrIter =
//        dictionaryConstructorTablePtr_->find(modelType);
//
//    if (cstrIter == dictionaryConstructorTablePtr_->end())
//    {
//        FatalErrorIn
//        (
//            "diffusivityModel::New(const dictionary&, "
//            "const volScalarField&, const concReactionThermo&, PtrList<label>&)"
//        )   << "Unknown diffusivityModel type "
//            << modelType << endl << endl
//            << "Valid diffusivityModels are : " << endl
//            << dictionaryConstructorTablePtr_->toc()
//            << exit(FatalError);
//    }
//
//    return autoPtr<diffusivityModel>
//        (cstrIter()(dict, thermo, Cstd));
//}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

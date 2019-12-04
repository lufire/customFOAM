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

//#include "diffusivityModel.H"
#include "electrolyteModel.H"
#include "conductivityModel.H"
//#include "volFields.H"

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(conductivityModel, 0);
defineRunTimeSelectionTable(conductivityModel, dictionary);
defineRunTimeSelectionTable(conductivityModel, constant);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
conductivityModel::conductivityModel
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    dict_(dict),
    mesh_(electrolyte.thermo().T().mesh()),
    thermo_(electrolyte.thermo()),
    species_(electrolyte.thermo().composition().species()),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "kappa",
            sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
            0.0
        )
    )
{}

conductivityModel::conductivityModel
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    dict_(dict),
    mesh_(thermo.T().mesh()),
    thermo_(thermo),
    species_(thermo.composition().species()),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "kappa",
            sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
            0.0
        )
    )
{}

//conductivityModel::conductivityModel
//(
//    const dictionary& dict,
//    const concReactionThermo& thermo,
//    const PtrList<label>& z
//)
//:
//    dict_(dict),
//    mesh_(thermo.T().mesh()),
//    thermo_(thermo),
//    species_(thermo.composition().species()),
//    kappa_
//    (
//        IOobject
//        (
//            "kappa",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar
//        (
//            "kappa",
//            sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
//            0.0
//        )
//    )
//{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<conductivityModel> conductivityModel::New
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
{
    word modelType(dict.lookup("typeName"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "conductivityModel::New(const dictionary&, "
            "const electrolyteModel&)"
        )   << "Unknown conductivityModel type "
            << modelType << endl << endl
            << "Valid conductivityModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<conductivityModel>
        (cstrIter()(dict, electrolyte));
}

autoPtr<conductivityModel> conductivityModel::New
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
{
    word modelType(dict.lookup("typeName"));

    typename constantConstructorTable::iterator cstrIter =
        constantConstructorTablePtr_->find(modelType);

    if (cstrIter == constantConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "conductivityModel::New(const dictionary&, "
            "const concReactionThermo&)"
        )   << "Unknown conductivityModel type "
            << modelType << endl << endl
            << "Valid conductivityModels are : " << endl
            << constantConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<conductivityModel>
        (cstrIter()(dict, thermo));
}

//autoPtr<conductivityModel> conductivityModel::New
//(
//    const dictionary& dict,
//    const concReactionThermo& thermo,
//    const PtrList<label>& z
//)
//{
//    word modelType(dict.lookup("typeName"));
//
//    typename dictionaryConstructorTable::iterator cstrIter =
//        dictionaryConstructorTablePtr_->find(modelType);
//
//    if (cstrIter == dictionaryConstructorTablePtr_->end())
//    {
//        FatalErrorIn
//        (
//            "conductivityModel::New(const dictionary&, "
//            "const volScalarField&, const concReactionThermo&, PtrList<label>&)"
//        )   << "Unknown conductivityModel type "
//            << modelType << endl << endl
//            << "Valid conductivityModels are : " << endl
//            << dictionaryConstructorTablePtr_->toc()
//            << exit(FatalError);
//    }
//
//    return autoPtr<conductivityModel>
//        (cstrIter()(dict, thermo, z));
//}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

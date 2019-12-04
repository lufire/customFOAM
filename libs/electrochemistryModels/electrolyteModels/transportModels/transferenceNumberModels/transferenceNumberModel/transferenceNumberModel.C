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
#include "transferenceNumberModel.H"
//#include "volFields.H"

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(transferenceNumberModel, 0);
defineRunTimeSelectionTable(transferenceNumberModel, dictionary);
defineRunTimeSelectionTable(transferenceNumberModel, constant);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
transferenceNumberModel::transferenceNumberModel
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    dict_(dict),
    mesh_(electrolyte.thermo().T().mesh()),
    species_(electrolyte.thermo().composition().species()),
    t_(species_.size())
{
    forAll(t_,i)
    {
        t_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "t_" + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("t", dimless, 0.0)
            )
        );
    }
}

transferenceNumberModel::transferenceNumberModel
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    dict_(dict),
    mesh_(thermo.T().mesh()),
    species_(thermo.composition().species()),
    t_(species_.size())
{
    forAll(t_,i)
    {
        t_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "t_" + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("t", dimless, 0.0)
            )
        );
    }
}

//transferenceNumberModel::transferenceNumberModel
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

autoPtr<transferenceNumberModel> transferenceNumberModel::New
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
            "transferenceNumberModel::New(const dictionary&, "
            "const electrolyteModel&)"
        )   << "Unknown transferenceNumberModel type "
            << modelType << endl << endl
            << "Valid transferenceNumberModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<transferenceNumberModel>
        (cstrIter()(dict, electrolyte));
}

autoPtr<transferenceNumberModel> transferenceNumberModel::New
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
            "transferenceNumberModel::New(const dictionary&, "
            "const concReactionThermo&)"
        )   << "Unknown transferenceNumberModel type "
            << modelType << endl << endl
            << "Valid transferenceNumberModels are : " << endl
            << constantConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<transferenceNumberModel>
        (cstrIter()(dict, thermo));
}

//autoPtr<transferenceNumberModel> transferenceNumberModel::New
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
//            "transferenceNumberModel::New(const dictionary&, "
//            "const volScalarField&, const concReactionThermo&, PtrList<label>&)"
//        )   << "Unknown transferenceNumberModel type "
//            << modelType << endl << endl
//            << "Valid transferenceNumberModels are : " << endl
//            << dictionaryConstructorTablePtr_->toc()
//            << exit(FatalError);
//    }
//
//    return autoPtr<transferenceNumberModel>
//        (cstrIter()(dict, thermo, z));
//}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

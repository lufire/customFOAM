/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "electrodeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

defineTypeNameAndDebug(electrodeModel, 0);
defineRunTimeSelectionTable(electrodeModel, electrolyte);
//defineRunTimeSelectionTable(electrodeModel, dictionary);
defineRunTimeSelectionTable(electrodeModel, myCfdemCloud);

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //


void electrodeModel::initializeMembers()
{
    word controlSwitch(electrodeDict_.lookup("currentControl"));
    if(controlSwitch == "true")
    {
        currentControl_ = true;
    }

    forAll(C_, speciesI)
    {
        stoichCoeff_.set
        (
            speciesI, new label 
            (
                readLabel
                (
                    electrodeDict_.subDict("stoichCoeffs").lookup(
                        electrolyte_.species()[speciesI])
                )
            )
        );
        reactionOrder_.set
        (
            speciesI, new label 
            (
                readLabel
                (
                    electrodeDict_.subDict("reactionOrders").lookup(
                        electrolyte_.species()[speciesI])
                )
            )
        );
        CRef_.set
        (
            speciesI, new scalar
            (
                readScalar
                (
                    electrodeDict_.subDict("CRefs").lookup(
                        electrolyte_.species()[speciesI])
                )
            )
        );
        if(electrolyte_.species()[speciesI] == redName_)
        {
            redID_ = speciesI;
        }
        else if(electrolyte_.species()[speciesI] == oxName_)
        {
            oxID_ = speciesI;
        }
    }
    if(oxID_ == -1 || redID_ == -1)
    {
        FatalErrorIn("electrodeModel::electrodeModel()")
            << "oxidant or reductant name does not"
            << " correspond to given species names"
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
electrodeModel::electrodeModel
(
    const dictionary& dict, 
    electrolyteModels::electrolyteModel& electrolyte 
)
:
    electrodeDict_(dict),
    electrolyte_(electrolyte),
    mesh_(electrolyte.thermo().T().mesh()),
    patchName_(dict.lookup("patchName")),
    patchID_(mesh_.boundary().findPatchID(patchName_)),
    currentControl_(false),
    surfaceDensity_(mesh_.V()*0.0),
    //surfaceDensity_
    //(
    //    IOobject
    //    (
    //        "surfaceDensity",
    //        mesh_.time().timeName(),
    //        mesh_,
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    mesh_,
    //    dimensionedScalar("surfaceDensity", dimArea/dimVolume, 0.0)
    //),
    T_(electrolyte.thermo().T()),
    oxName_(electrodeDict_.lookup("oxidant")),
    redName_(electrodeDict_.lookup("reductant")),
    C_(electrolyte.thermo().composition().C()),
    electronNumber_(readLabel(electrodeDict_.lookup("electronNumber"))),
    stoichCoeff_(C_.size()),
    reactionOrder_(C_.size()),
    CRef_(C_.size()),
    iEx_(readScalar(electrodeDict_.lookup("exchangeCurrentDensity"))),
    eqPotential_(readScalar(
        electrodeDict_.lookup("equilibriumPotential"))),
    alphaA_(readScalar(electrodeDict_.lookup("alphaA"))),
    alphaC_(readScalar(electrodeDict_.lookup("alphaC"))),
    iTrans_
    (
        IOobject
        (
            "iTrans",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("iTrans", dimCurrent/dimVolume, 0.0)
    ),
    is_
    (
        IOobject
        (
            "is",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("is", dimCurrent/dimArea, vector(0.0, 0.0, 0.0))
    ),
    oxID_(-1),
    redID_(-1),
    sigma_
    (
        "sigma",
        sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
        readScalar(electrodeDict_.lookup("sigma"))
    ),
    sigmaEff_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigma_
    )
{
    initializeMembers();
}

//electrodeModel::electrodeModel
//(
//    const dictionary& electrodeDict, 
//    const concReactionThermo& thermo
//)
//:
//    electrodeDict_(electrodeDict),
//    mesh_(thermo.T().mesh()),
//    patchName_(electrodeDict.lookup("patchName")),
//    patchID_(mesh_.boundary().findPatchID(patchName_)),
//    surfaceDensity_(mesh_.V()*0.0),
//    //surfaceDensity_
//    //(
//    //    IOobject
//    //    (
//    //        "surfaceDensity",
//    //        mesh_.time().timeName(),
//    //        mesh_,
//    //        IOobject::NO_READ,
//    //        IOobject::AUTO_WRITE
//    //    ),
//    //    mesh_,
//    //    dimensionedScalar("surfaceDensity", dimArea/dimVolume, 0.0)
//    //),
//    T_(thermo.T()),
//    oxName_(electrodeDict_.lookup("oxidant")),
//    redName_(electrodeDict_.lookup("reductant")),
//    C_(thermo.composition().C()),
//    electronNumber_(readLabel(electrodeDict_.lookup("electronNumber"))),
//    stoichCoeff_(C_.size()),
//    reactionOrder_(C_.size()),
//    CRef_(C_.size()),
//    iEx_(readScalar(electrodeDict_.lookup("exchangeCurrentDensity"))),
//    eqPotential_(readScalar(
//        electrodeDict_.lookup("equilibriumPotential"))),
//    alphaA_(readScalar(electrodeDict_.lookup("alphaA"))),
//    alphaC_(readScalar(electrodeDict_.lookup("alphaC"))),
//    iTrans_
//    (
//        IOobject
//        (
//            "iTrans",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar("iTrans", dimCurrent/dimVolume, 0.0)
//    ),
//    is_
//    (
//        IOobject
//        (
//            "is",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh_,
//        dimensionedVector("is", dimCurrent/dimArea, vector(0.0, 0.0, 0.0))
//    ),
//    oxID_(-1),
//    redID_(-1),
//    sigma_
//    (
//        "sigma",
//        sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
//        readScalar(electrodeDict_.lookup("sigma"))
//    ),
//    sigmaEff_
//    (
//        IOobject
//        (
//            "sigma",
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh_,
//        sigma_
//    )
//{
//    initializeMembers();
//}


electrodeModel::electrodeModel
(
    const dictionary& dict, 
    electrolyteModels::electrolyteModel& electrolyte,
    const myCfdemCloud& particleCloud
)
:
    electrodeDict_(dict),
    electrolyte_(electrolyte),
    mesh_(electrolyte.thermo().T().mesh()),
    patchName_(dict.lookup("patchName")),
    patchID_(mesh_.boundary().findPatchID(patchName_)),
    currentControl_(false),
    surfaceDensity_(mesh_.V()*0.0),
    //surfaceDensity_
    //(
    //    IOobject
    //    (
    //        "surfaceDensity",
    //        mesh_.time().timeName(),
    //        mesh_,
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    mesh_,
    //    dimensionedScalar("surfaceDensity", dimArea/dimVolume, 0.0)
    //),
    T_(electrolyte.thermo().T()),
    oxName_(electrodeDict_.lookup("oxidant")),
    redName_(electrodeDict_.lookup("reductant")),
    C_(electrolyte.thermo().composition().C()),
    electronNumber_(readLabel(electrodeDict_.lookup("electronNumber"))),
    stoichCoeff_(C_.size()),
    reactionOrder_(C_.size()),
    CRef_(C_.size()),
    iEx_(readScalar(electrodeDict_.lookup("exchangeCurrentDensity"))),
    eqPotential_(readScalar(
        electrodeDict_.lookup("equilibriumPotential"))),
    alphaA_(readScalar(electrodeDict_.lookup("alphaA"))),
    alphaC_(readScalar(electrodeDict_.lookup("alphaC"))),
    iTrans_
    (
        IOobject
        (
            "iTrans",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("iTrans", dimCurrent/dimVolume, 0.0)
    ),
    is_
    (
        IOobject
        (
            "is",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("is", dimCurrent/dimArea, vector(0.0, 0.0, 0.0))
    ),
    oxID_(-1),
    redID_(-1),
    sigma_
    (
        "sigma",
        sqr(dimCurrent)*pow3(dimTime)/(dimVol*dimMass),
        readScalar(electrodeDict_.lookup("sigma"))
    ),
    sigmaEff_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigma_
    )
{
    initializeMembers();
}
// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

electrodeModel::~electrodeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam
// ************************************************************************* //

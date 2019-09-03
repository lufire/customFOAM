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

#include "electrochemistryModel.H"
//#include "diffusivityModel.H"
//#include "conductivityModel.H"
//#include "electrodeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{

defineTypeNameAndDebug(electrochemistryModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

electrochemistryModel::electrochemistryModel
(
    const concReactionThermo& thermo,
    const myCfdemCloud& particleCloud
)
:
    mesh_(thermo.T().mesh()),
    electrochemDict_ 
    (
        IOobject
        (
            "electrochemicalProperties",
            mesh_.time().caseConstant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    isActive_(false),
    currentControl_(false),
    electrolyte_
    (
        electrolyteModels::electrolyteModel::New
        (
            electrochemDict_.subDict("electrolyteModel"),
            thermo
        )
    ),
    electrodeTable_(electrochemDict_.lookup("electrodes")),
    electrodes_(electrodeTable_.size())
{
    word electrochemSwitch(electrochemDict_.lookup("electrochemistry"));
    if(electrochemSwitch == "on")
    {
        isActive_ = true;
    }
    word controlSwitch(electrochemDict_.lookup("currentControl"));
    if(controlSwitch == "on")
    {
        isActive_ = true;
    }
    //Create electrode models
    forAll(electrodeTable_, electrodeI)
    {
        dictionary electrodeDict
        (
            electrochemDict_.subDict(electrodeTable_[electrodeI])
        );
        const word modelType(electrodeDict.lookup("typeName"));
        Info << "electrode typeName " << modelType << endl;
        if (modelType == "cfdemPorousElectrode")
        {
            electrodes_.set
            (
                electrodeI, 
                electrodeModels::electrodeModel::New
                (
                    electrodeDict, 
                    electrolyte_(),
                    particleCloud
                )
            );
        }
        else
        {
            electrodes_.set
            (
                electrodeI, 
                electrodeModels::electrodeModel::New
                (
                    electrodeDict, 
                    electrolyte_() 
                )
            );
        }
    }
}

electrochemistryModel::electrochemistryModel
(
    const concReactionThermo& thermo
)
:
    mesh_(thermo.T().mesh()),
    electrochemDict_ 
    (
        IOobject
        (
            "electrochemicalProperties",
            mesh_.time().caseConstant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    isActive_(false),
    currentControl_(false),
    electrolyte_
    (
        electrolyteModels::electrolyteModel::New
        (
            electrochemDict_.subDict("electrolyteModel"),
            thermo
        )
    ),
    electrodeTable_(electrochemDict_.lookup("electrodes")),
    electrodes_(electrodeTable_.size())
{
    word electrochemSwitch(electrochemDict_.lookup("electrochemistry"));
    if(electrochemSwitch == "on")
    {
        isActive_ = true;
    }
    word controlSwitch(electrochemDict_.lookup("currentControl"));
    if(controlSwitch == "on")
    {
        isActive_ = true;
    }
    //Create electrode models
    forAll(electrodeTable_, electrodeI)
    {
        dictionary electrodeDict
        (
            electrochemDict_.subDict(electrodeTable_[electrodeI])
        );
        electrodes_.set
        (
            electrodeI, 
            electrodeModels::electrodeModel::New
            (
                electrodeDict, 
                electrolyte_() 
            )
        );
    }
}
// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

electrochemistryModel::~electrochemistryModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void electrochemistryModel::correct()
{
    electrolyte_->correct();
    forAll(electrodes_,electrodeI)
    {
        electrodes_[electrodeI].correct();
    }
}

// ************************************************************************* //

} // End namespace electrochemistryModels 
} // End namespace Foam

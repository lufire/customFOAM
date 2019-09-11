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

#include "NernstPlanck.H"
#include "addToRunTimeSelectionTable.H"
#include "electrolyteModel.H"
#include "fvCFD.H"
//#include "fvm.H"
//#include "diffusivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(NernstPlanck, 0);
addToRunTimeSelectionTable
(
    conductivityModel,
    NernstPlanck, 
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NernstPlanck::NernstPlanck
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    conductivityModel(dict, electrolyte),
    balanceIndex_(0),
    z_(electrolyte.z()),
    C_(electrolyte.thermo().composition().C()),
    D_(electrolyte.D())
{
    update();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void NernstPlanck::update()
{
    const dimensionedScalar& R = constant::physicoChemical::R;
    const dimensionedScalar& F = constant::physicoChemical::F;
    const volScalarField& T = thermo_.T();
    kappa_ *= 0.0;
    forAll(C_, speciesI)
    {
        kappa_ += F*F/(R*T)*z_[speciesI]*z_[speciesI]*D_[speciesI]*C_[speciesI];
        //kappa_ *= F*F/(R*T);
    }
}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

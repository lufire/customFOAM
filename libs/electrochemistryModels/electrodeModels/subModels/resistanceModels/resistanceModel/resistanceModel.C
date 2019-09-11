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

#include "resistanceModel.H"
//#include "conductivityModel.H"
#include "electrodeModel.H"
//#include "volFields.H"

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

defineTypeNameAndDebug(resistanceModel, 0);
defineRunTimeSelectionTable(resistanceModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
resistanceModel::resistanceModel
(
    const dictionary& dict,
    const electrodeModel& electrode 
)
:
    dict_(dict),
    electrode_(electrode),
    thermo_(electrode.electrolyte().thermo()),
    mesh_(thermo_.T().mesh()),
    species_(thermo_.composition().species()),
    Cs_(species_.size()),
    phiEls_
    (
        electrode.electrolyte().phiE().boundaryField()[electrode.patchIndex()]
    )
    //rotationMatrix_(mesh_.boundary()[electrode.patchName()].nf()->size())
{
    forAll(species_, specieI)
    {
        const volScalarField& Ci = thermo_.composition().C()[specieI];
        Cs_.set
        (
            specieI, new scalarField(Ci.boundaryField()[electrode.patchIndex()])
        );
    }

    //const label patchName = electrode_.patchName();
    //const vectorField& n = mesh_.boundary()[patchName].nf()();
    //const scalarField theta = acos(vector(1.0,0,0)&n);
    //const vectorField axis = vector(1.0,0,0)^n;
    //const scalarField a = cos(theta*0.5);
    //const scalarField b = (axis&vector(1.0,0,0))*sin(theta*0.5);
    //const scalarField c = (axis&vector(0,1.0,0))*sin(theta*0.5);
    //const scalarField d = (axis&vector(0,0,1.0))*sin(theta*0.5);
    //forAll(rotationMatrix_, i)
    //{
    //    if (theta[i] == 0.0)
    //    {    
    //        rotationMatrix_[i] 
    //            = tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0); 
    //    }
    //    else if (theta[i] == constant::mathematical::pi)
    //    {
    //        rotationMatrix_[i] 
    //            = tensor(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0); 
    //    }
    //    else
    //    {
    //        rotationMatrix_[i] = tensor
    //        (
    //            a[i]*a[i]+b[i]*b[i]-c[i]*c[i]-d[i]*d[i], 
    //            2.0*(b[i]*c[i]-a[i]*d[i]), 
    //            2.0*(b[i]*d[i]+a[i]*c[i]),
    //            2.0*(b[i]*c[i]+a[i]*d[i]), 
    //            a[i]*a[i]+c[i]*c[i]-b[i]*b[i]-d[i]*d[i], 
    //            2.0*(c[i]*d[i]-a[i]*b[i]),
    //            2.0*(b[i]*d[i]-a[i]*c[i]), 
    //            2.0*(c[i]*d[i]+a[i]*b[i]), 
    //            a[i]*a[i]+d[i]*d[i]-b[i]*b[i]-c[i]*d[i]
    //        );
    //    }
    //}

}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<resistanceModel> resistanceModel::New
(
    const dictionary& dict,
    const electrodeModel& electrode 
)
{
    word modelType(dict.lookup("typeName"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "resistanceModel::New(const dictionary&, "
            "const electrodeModel&)"
        )   << "Unknown resistanceModel type "
            << modelType << endl << endl
            << "Valid resistanceModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<resistanceModel>
        (cstrIter()(dict, electrode));
}

} // End namespace electrodeModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

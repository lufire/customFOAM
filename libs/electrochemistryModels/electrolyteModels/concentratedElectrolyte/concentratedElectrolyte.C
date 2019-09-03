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

#include "concentratedElectrolyte.H"
#include "diffusivityModel.H"
#include "conductivityModel.H"

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels 
{
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
concentratedElectrolyte<ThermoType>::
concentratedElectrolyte
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    electrolyteModel(dict, thermo),
    speciesThermo_
    (
        dynamic_cast<const concMultiComponentMixture<ThermoType>&>
            (this->thermo_.composition()).speciesData()
    )
{   
    //D_.setSize(species().size()-1);
    
    updateMolarFractions();
    //updateCoefficients();    
    //forAll(D_, i)
    //{
    //    D_.set
    //    (
    //        i, new volScalarField
    //        (
    //            IOobject
    //            (
    //                "D_" + species()[i],
    //                mesh_.time().timeName(),
    //                mesh_,
    //                IOobject::NO_READ,
    //                IOobject::NO_WRITE
    //            ),
    //            mesh_,
    //            dimensionedScalar("D", dimensionSet(0, 2, -1, 0, 0), 0.0)
    //        )
    //    );
    //} 

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void concentratedElectrolyte<ThermoType>::
updateCoefficients()
{     
    diffModel_->update();
    kappaModel_->update();
} 


template<class ThermoType>
void concentratedElectrolyte<ThermoType>::correct()
{     
    updateMolarFractions();
    updateCoefficients();
} 

//template<class ThermoType>
//Foam::scalar concentratedElectrolyte<ThermoType>::correct
//(
//    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
//)
//{
//    updateCoefficients();
//
//    scalar maxResidual = 0;
//    scalar eqnResidual = 1;
//
//    volScalarField yt = 0.0*thermo_.composition().Y(0);
//    surfaceScalarField nt = turbulence_.phi();
//    
//    forAll(this->D_, i)
//    {  
//        volScalarField& yi = thermo_.composition().Y(i);
//        surfaceScalarField& ni = n_[i];
//
//        tmp<fv::convectionScheme<scalar> > mvConvection
//        (
//            fv::convectionScheme<scalar>::New
//            (
//                mesh_,
//                fields,
//                turbulence_.phi(),
//                mesh_.divScheme("div(phi,Yi_h)")
//            )
//        );
//
//        if (mesh_.relaxField("Yi"))//Mohsen
//        {
//            yi.storePrevIter();
//        }
//            
//        tmp<fvScalarMatrix> yEqn
//        (   
//            fvm::ddt(thermo_.rho(), yi)
////           + fvm::div(turbulence_.phi(), yi, "div(phi,Yi_h)")
//          + mvConvection->fvmDiv(turbulence_.phi(), yi)
//          - fvm::laplacian(D_[i],yi, "laplacian(D,Yi)")
//          ==
//            Sy_[i]
//        );
//
//        eqnResidual = solve(yEqn() , mesh_.solver("Yi")).initialResidual();
//        maxResidual = max(eqnResidual, maxResidual);
//
//        if (mesh_.relaxField("Yi"))//Mohsen
//        {
//	  yi.relax(mesh_.fieldRelaxationFactor("Yi"));//Mohsen
//        }
//
//        yi.max(0.0);
////         yi.min(1.0);
//
//        ni = yEqn().flux();
//
//        nt -= ni;
//        yt += yi;  
//    }
//     
//    // Calculate inert species
//    volScalarField& yInert = thermo_.composition().Y()[inertIndex_];
//    yInert = 1 - yt;
//    forAll(yInert.boundaryField(), boundaryI)
//    {
//        forAll(yInert.boundaryField()[boundaryI], faceI)
//        {
//            yInert.boundaryField()[boundaryI][faceI] = 1- yt.boundaryField()[boundaryI][faceI];
//        }
//    }
//    yInert.max(0.0);
//    n_[inertIndex_] = nt;
//          
//    updateMolarFractions();
//
//    return maxResidual;
//}

//template<class ThermoType>
//bool concentratedElectrolyte<ThermoType>::read()
//{
//    if (regIOobject::read())
//    {
//        return true;
//    }
//    else
//    {
//        return false;
//    }
//}
   
} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

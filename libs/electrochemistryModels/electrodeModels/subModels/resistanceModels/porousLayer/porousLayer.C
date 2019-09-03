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

#include "porousLayer.H"
#include "electrodeModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels 
{

defineTypeNameAndDebug(porousLayer, 0);
addToRunTimeSelectionTable
(
    resistanceModel,
    porousLayer, 
    dictionary 
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

porousLayer::porousLayer
(
    const dictionary& dict,
    const electrodeModel& electrode 
)
:
    resistanceModel(dict, electrode),
    eps_(readScalar(dict.lookup("initialPorosity"))),
    delta_(readScalar(dict.lookup("initialThickness"))),
    deltaField_
    (
        (mesh_.boundary()[electrode.patchName()].magSf()
       /mesh_.boundary()[electrode.patchName()].magSf())*delta_
    ),
    brugg_(readScalar(dict.lookup("BruggemannCoeff")))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void porousLayer::correctPotential()
{
    const label patchID = electrode_.patchIndex();
    const word patchName = electrode_.patchName();

    const volVectorField& i =
        electrode_.electrolyte().i();

    const volVectorField& iD =
        electrode_.electrolyte().iD();


    const volScalarField& phiEl = electrode_.electrolyte().phiE();
    const volScalarField& kappa =
        electrode_.electrolyte().kappaModel().kappa();

    // Transform vectors into local coordinate system 
    //vectorField NoxLocal = Noxp; 
    //vectorField NredLocal = Nredp; 
    //vectorField iLocal = ip; 

    //forAll(rotationMatrix_, i)
    //{
    //    NoxLocal[i] = rotationMatrix_[i]&Noxp[i];
    //    NredLocal[i] = rotationMatrix_[i]&Nredp[i];
    //    iLocal[i] = rotationMatrix_[i]&ip[i];
    //}
    // 
    //const scalarField NoxNormal = NoxLocal&vector(1.0,0,0);
    //const scalarField NredNormal = NredLocal&vector(1.0,0,0);
    //const scalarField iNormal = iLocal&vector(1.0,0,0);

    if (i.dimensions() == dimCurrent/dimArea)
    {
        if(delta_ > 0.0)
        {
            surfaceScalarField snGrad = 
                (
                    -(fvc::interpolate(i) & mesh_.Sf())
                    +(fvc::interpolate(iD) & mesh_.Sf())
                )/(mesh_.magSf()*fvc::interpolate(kappa,"interpolate(kappa)"));

            const scalarField& phiElp = phiEl.boundaryField()[patchID];
            const scalarField& snGradp = snGrad.boundaryField()[patchID];
            phiEls_ = phiElp + snGradp*deltaField_;
            //Info << "phiEls_[0] = " << phiEls_[0] << endl;
            //Info << "phiElp[0] = " << phiElp[0] << endl;

        }
        else
        {
            const scalarField& phiElp = phiEl.boundaryField()[patchID];
            phiEls_ = phiElp ;
        }

    }
    else
    {
        FatalErrorIn(
                "resistanceModel::porousLayer::update()")
            << "dimensions of field "
            << i.name() 
            << "are not correct"
            << exit(FatalError);
    }
}

void porousLayer::correctSpecies()
{
    const label patchID = electrode_.patchIndex();
    const word patchName = electrode_.patchName();

    const PtrList<volScalarField>& C = 
        thermo_.composition().C();

    const PtrList<volVectorField>& N = 
        electrode_.electrolyte().N();

    const volVectorField& i =
        electrode_.electrolyte().i();

    const PtrList<volScalarField>& t =
        electrode_.electrolyte().kappaModel().t();

    const PtrList<label>& z =
        electrode_.electrolyte().z();

    const PtrList<volScalarField>& D =
        electrode_.electrolyte().diffModel().D();

    // Transform vectors into local coordinate system 
    //vectorField NoxLocal = Noxp; 
    //vectorField NredLocal = Nredp; 
    //vectorField iLocal = ip; 

    //forAll(rotationMatrix_, i)
    //{
    //    NoxLocal[i] = rotationMatrix_[i]&Noxp[i];
    //    NredLocal[i] = rotationMatrix_[i]&Nredp[i];
    //    iLocal[i] = rotationMatrix_[i]&ip[i];
    //}
    // 
    //const scalarField NoxNormal = NoxLocal&vector(1.0,0,0);
    //const scalarField NredNormal = NredLocal&vector(1.0,0,0);
    //const scalarField iNormal = iLocal&vector(1.0,0,0);

    if (N[0].dimensions() == dimMoles/dimArea/dimTime)
    {
        if(delta_ > 0.0)
        {
            const dimensionedScalar& F = constant::physicoChemical::F;
            forAll(C,specieI)
            {
                surfaceScalarField snGrad =
                    (
                        -(mesh_.Sf() & fvc::interpolate(N[specieI]))
                    );
                if(z[specieI] != 0)
                {
                    snGrad +=
                        (
                            -mesh_.Sf() 
                             & (fvc::interpolate
                                (
                                    t[specieI]/(z[specieI]*F),
                                    "interpolate(t)"
                                )
                             *fvc::interpolate(i))
                        );
       
                }

                snGrad /=
                    mesh_.magSf()*fvc::interpolate
                    (
                        D[specieI]*pow(eps_,brugg_),
                        "interpolate(D)"
                    );

                const scalarField& Cp = C[specieI].boundaryField()[patchID];
                const scalarField& snGradp = snGrad.boundaryField()[patchID];
                
                Cs_[specieI] = Cp + snGradp*deltaField_;

                forAll(Cs_[specieI],faceI)
                {
                    if(Cs_[specieI][faceI]<0.0)
                    {
                        Cs_[specieI][faceI] = 0.0;
                    }
                }
                //Info << "Cs_[" << specieI <<"][5] = " << Cs_[specieI][5] << endl;
                //Info << "Cp[" << specieI << "][5] = " << Cp[5] << endl;
            }

        }
        else
        {
            forAll(C,specieI)
            {
                const scalarField& Cp = C[specieI].boundaryField()[patchID];

                Cs_[specieI] = Cp;
            }
        }

    }
    else
    {
        FatalErrorIn(
                "resistanceModel::porousLayer::update()")
            << "dimensions of field "
            << N[0].name() 
            << "are not correct"
            << exit(FatalError);
    }
}

void porousLayer::correctSpeciesFlux(vectorField& N)
{
    const label patchID = electrode_.patchIndex();
    const word patchName = electrode_.patchName();
    const tmp<vectorField> tn = mesh_.boundary()[patchName].nf();
    const vectorField& n = tn(); 

    // Hard coded precipitation reaction
    const PtrList<volScalarField>& C = 
        thermo_.composition().C();
    const scalarField& Cref = C[2].boundaryField()[patchID];

    const PtrList<volScalarField>& D =
        electrode_.electrolyte().diffModel().D();
    const scalarField& Dp = D[0].boundaryField()[patchID];

    scalar Cstd = electrode_.electrolyte().Cstd().value();
    scalarField Csat = Cref;

    forAll(Csat,faceI)
    {
        if(Cref[faceI]> 2.1*Cstd)
        {
            Csat[faceI] = -0.21*Cstd + 0.975e-1*Cref[faceI]
                + 0.125e-2*(sqr(Cref[faceI])/Cstd);
        }
        else
        {
            Csat[faceI] = 0.0;
        }
    }
        
    N += Dp*pow(eps_,brugg_)*(Cs_[0]-Csat)/deltaField_*n*-1.0;
}

} // End namespace electrodeModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

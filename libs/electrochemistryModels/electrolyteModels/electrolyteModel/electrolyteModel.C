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

#include "electrolyteModel.H"
#include "dimensionedConstants.H"
#include "constants.H" 
#include "concChemistryModel.H"
//#include "diffusivityModel.H"
//#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels 
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(electrolyteModel, 0);
defineRunTimeSelectionTable(electrolyteModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

electrolyteModel::electrolyteModel
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    electrolyteDict_(dict),
    mesh_(thermo.T().mesh()),
    pZones_
    (
        mesh_,
        IOdictionary
        (
            IOobject
            (
                "porousZones",
                mesh_.time().constant(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    thermo_(thermo),
    species_(thermo.composition().species()),
    balanceSpecies_(dict.lookup("balanceSpecies")),
    inertSpecies_(dict.lookup("inertSpecies")),
    Cstd_(dict.lookup("standardConcentration")),
    epsilon_(dict.lookup("relativePermittivity")),
    x_(species_.size()),
    N_(species_.size()),
    S_(species_.size()),
    z_(species_.size()),
    diffModel_(),
    kappaModel_(),
    transNumberModel_(),
    phiE_
    (
        IOobject
        (
            "phiE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    i_
    (
        IOobject
        (
            "i",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "i", dimCurrent/dimArea, vector(0.0, 0.0, 0.0)
        )
    ),
    iD_
    (
        IOobject
        (
            "iD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "iD", dimCurrent/dimArea, vector(0.0, 0.0, 0.0)
        )
    )
        
{
    forAll(species_, i)
    {
        x_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "x_" + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimless
            )
        );

        N_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "N_" + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "N", dimMoles/(dimArea*dimTime), vector(0.0, 0.0, 0.0)
                )
            )
        );
        S_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "S_" + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("S", dimMoles/dimTime/dimVolume, 0.0)
            )
        );

        z_.set(i, new label(0));        
        const word& speciesName = species_[i];
        size_t found = speciesName.find("-");
        if(found == std::string::npos)
        {
            found = speciesName.find("+");
            if(found == std::string::npos)
            {
                z_[i] = 0;
            }
            else
            {
                z_[i] = speciesName.size() - found;
            }
        }
        else
        {
            z_[i] = (speciesName.size() - found)*(-1);
        }
        Info << "Name of Species: " << speciesName << ", z: " << z_[i] << endl;
    }

    // Construct diffusivity and conducitivty Models
    diffModel_.set
    (
        diffusivityModel::New
        (
            dict.subDict("diffusivityModel"),
            *this
        ).ptr()
    );
    kappaModel_.set
    (
        conductivityModel::New
        (
            dict.subDict("conductivityModel"),
            *this
        ).ptr()
    );
    transNumberModel_.set
    (
        transferenceNumberModel::New
        (
            dict.subDict("transferenceNumberModel"),
            *this
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//Foam::surfaceScalarField electrolyteModel::j
//(
//    const Foam::label i
//)
//{
//    return n(i) -
//        fvc::interpolate(thermo_.composition().Y(i)) * turbulence_.phi();
//}

void electrolyteModel::updateMolarFractions()
//void electrolyteModel::updateMolarFractions()
{
    const PtrList<volScalarField>& Y = thermo_.composition().Y(); 
    forAll(Y[0].internalField(), cellI)
    {
	    scalar tW = 0.0;
        forAll(Y, i)
        {
	        tW +=  Y[i].internalField()[cellI]/W(i);
        }
        tW = 1.0/tW;
        forAll(Y, i)
        {
	        x_[i].internalField()[cellI] =  Y[i].internalField()[cellI]*tW/W(i);
        }
    }

    forAll(Y[0].boundaryField(), patchI)
    {

	    forAll(Y[0].boundaryField()[patchI], faceI)
	    {
	        scalar ptW = 0.0;
		
            forAll(Y, i)
            {
	            ptW +=  Y[i].boundaryField()[patchI][faceI]/W(i);
	        }
	        ptW = 1.0/ptW;
            forAll(Y, i)
            {
	            x_[i].boundaryField()[patchI][faceI] =  Y[i].boundaryField()[patchI][faceI]*ptW/W(i);
	        }
        }
    }
}

//Foam::scalar electrolyteModel::correct
//(
//    const PtrList<volScalarField>& Y,
//    const volScalarField& gamma,
//    const concChemistryModel& chemistry,
//    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
//)
//{
//
//    forAll(Sy_, i)
//    {
//          Sy_[i] = gamma*RR(Y[i].name(), chemistry,i);
//    }
//
//    return correct(fields);
//}
//
//Foam::tmp<Foam::volScalarField> electrolyteModel::RR
//(
//    const word& Yname,
//    const concChemistryModel& chemistry,
//    const label& i
//) const
//{
//
//    tmp<volScalarField> tRR
//    (
//        new volScalarField
//        (
//            IOobject
//            (
//                "RR(" + Yname + ')',
//                chemistry.time().timeName(),
//                chemistry.mesh(),
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            chemistry.mesh(),
//            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
//            zeroGradientFvPatchScalarField::typeName
//        )
//    );
//
//    if (chemistry.chemistry())
//    {
//        tRR().internalField() = chemistry.RR(i);
//        tRR().correctBoundaryConditions();
//    }
//    return tRR;
//}

//bool electrolyteModel::read()
//{
//    return regIOobject::read();
//}

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

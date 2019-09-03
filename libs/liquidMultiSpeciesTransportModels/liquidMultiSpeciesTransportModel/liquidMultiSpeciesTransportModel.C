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

#include "liquidMultiSpeciesTransportModel.H"
#include "dimensionedConstants.H"
#include "constants.H" 
#include "concChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(liquidMultiSpeciesTransportModel, 0);
    defineRunTimeSelectionTable(liquidMultiSpeciesTransportModel, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidMultiSpeciesTransportModel::liquidMultiSpeciesTransportModel
(
    concReactionThermo& thermo,
    const incompressible::turbulenceModel& turbulence
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    mesh_(thermo.T().mesh()),

    // M.Lindner, A.Alexiou 16.12.2014
    pZones_
    (
        thermo.T().mesh(),
        IOdictionary
        (
            IOobject
            (
                "porousZones",
                thermo.T().mesh().time().constant(),
                thermo.T().mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),

    thermo_(thermo),

    // The last specie in the list is the default inertspecie
    inertIndex_(species().size()-1),

    turbulence_(turbulence)
{
    // Construct the diffusivity model
    DijModel_.set(new diffusivityModel(thermo.p(), thermo.T(), pZones_, species()));

    x_.setSize(species().size());
    n_.setSize(species().size());
    Sy_.setSize(species().size());
    z_.setSize(species().size()); 
    
    forAll(thermo.composition().Y(), i)
    {
        x_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "x_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimless
            )
        );

        n_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "n_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("n", dimMass/dimTime, 0.0)
            )
        );

        Sy_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "Sy_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("Sy", dimMass/dimTime/dimVolume, 0.0)
            )
        );
        
        z_.set(i, new label(0));        
        
        const word& speciesName = species()[i];
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::surfaceScalarField Foam::liquidMultiSpeciesTransportModel::j
(
    const Foam::label i
)
{
    return n(i) -
        fvc::interpolate(thermo_.composition().Y(i)) * turbulence_.phi();
}

void Foam::liquidMultiSpeciesTransportModel::updateMolarFractions()
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

Foam::scalar Foam::liquidMultiSpeciesTransportModel::correct
(
    const PtrList<volScalarField>& Y,
    const volScalarField& kappa,
    const concChemistryModel& chemistry,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{

    forAll(Sy_, i)
    {
          Sy_[i] = kappa*RR(Y[i].name(), chemistry,i);
    }

    return correct(fields);
}

Foam::tmp<Foam::volScalarField> Foam::liquidMultiSpeciesTransportModel::RR
(
    const word& Yname,
    const concChemistryModel& chemistry,
    const label& i
) const
{

    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + Yname + ')',
                chemistry.time().timeName(),
                chemistry.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            chemistry.mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (chemistry.chemistry())
    {
        tRR().internalField() = chemistry.RR(i);
        tRR().correctBoundaryConditions();
    }
    return tRR;
}

bool Foam::liquidMultiSpeciesTransportModel::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more dnuils.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PotaninStructure.H"
#include "addToRunTimeSelectionTable.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PotaninStructure, 0);
    addToRunTimeSelectionTable(structureModel, PotaninStructure, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PotaninStructure::PotaninStructure
(
    const word& name,
    const dictionary& dict,
    const rheologyModel& rheology,
    const volVectorField& U
)
:
    structureModel(name, rheology, U),
    dict_(dict),
    tR_(dict.lookup("tR")),
    tB_(dict.lookup("tB")),
    srRef_(dict.lookup("srRef")),
    nR_(dict.lookup("nR")),
    nB_(dict.lookup("nB")),
    mR_(dict.lookup("mR")),
    mB_(dict.lookup("mB")),
    f_
    (
        IOobject
        (
            "f",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        rheology.nuInf()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PotaninStructure::correct(const Foam::surfaceScalarField& phi)
{
    const volScalarField lambdaEq = rheology_.nuEq()/f_;
    //const volScalarField lambdaEq = rheology_.nuEq()/rheology_.nuEq();
    lambda_.max(SMALL);
    
    const volScalarField& sr = rheology_.sr();

    expLambdaSource_ = 0.0*lambda_/tR_;
    volScalarField n = lambda_*0.0;

    dimensionedScalar tc = tR_;
    forAll(expLambdaSource_, cellI)
    {
        if(lambda_[cellI] < lambdaEq[cellI])
        {
            tc = tR_*pow(srRef_.value()/sr[cellI], mR_.value());
            expLambdaSource_[cellI] = 1.0/tc.value();
            n[cellI] = nR_.value();
        }
        else
        {
            tc = tB_*pow(srRef_.value()/sr[cellI], mB_.value());
            expLambdaSource_[cellI] = -1.0/tc.value();
            n[cellI] = nB_.value();
        }
    }
    forAll(expLambdaSource_.boundaryField(), patchI)
    {
        forAll(expLambdaSource_.boundaryField()[patchI], faceI)
        {
            if(lambda_.boundaryField()[patchI][faceI] 
                < lambdaEq.boundaryField()[patchI][faceI])
            {
                tc = tR_*pow(srRef_.value()/sr.boundaryField()[patchI][faceI], 
                        mR_.value());
                expLambdaSource_.boundaryField()[patchI][faceI] = 1.0/tc.value();
                n.boundaryField()[patchI][faceI] = nR_.value();
            }
            else
            {
                tc = tB_*pow(srRef_.value()/sr.boundaryField()[patchI][faceI],
                        mB_.value());
                expLambdaSource_.boundaryField()[patchI][faceI] = -1.0/tc.value();
                n.boundaryField()[patchI][faceI] = nB_.value();
            }
        }
    }
    expLambdaSource_ *= pow(mag(lambda_ - lambdaEq), n);

    structureModel::correct(phi);

}


// ************************************************************************* //

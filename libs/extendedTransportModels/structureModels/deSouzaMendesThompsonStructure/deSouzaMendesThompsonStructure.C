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

#include "deSouzaMendesThompsonStructure.H"
#include "addToRunTimeSelectionTable.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deSouzaMendesThompsonStructure, 0);
    addToRunTimeSelectionTable(structureModel, deSouzaMendesThompsonStructure, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deSouzaMendesThompsonStructure::deSouzaMendesThompsonStructure
(
    const word& name,
    const dictionary& dict,
    const rheologyModel& rheology,
    const volVectorField& U
)
:
    structureModel(name, rheology, U),
    dict_(dict),
    teq_(dict.lookup("teq")),
    a_(dict.lookup("a")),
    b_(dict.lookup("b"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::deSouzaMendesThompsonStructure::correct(const Foam::surfaceScalarField& phi)
{
    const volScalarField lambdaEq = log(rheology_.nuEq()/rheology_.nuInf());
    const dimensionedScalar lambda0 = 
        min
        (
            log(rheology_.nu0()/rheology_.nuInf()), 
            dimensionedScalar("lambdaMax", dimless, 1e6)
        );
    lambda_.max(SMALL);

    expLambdaSource_ = 
         (pow((1.0/lambda_-1.0/lambda0), a_)
         - pow((lambda_/lambdaEq),b_)*pow((1.0/lambdaEq - 1.0/lambda0), a_))
         /(teq_);
    //dimensionedScalar diffStab("diffStab", dimArea/(dimTime), 0.5);

    structureModel::correct(phi);
}


// ************************************************************************* //

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
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "structureModel.H"
#include "solutionControl.H"
#include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structureModel, 0);
    defineRunTimeSelectionTable(structureModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureModel::structureModel
(
    const word& name,
    const rheologyModel& rheology,
    const volVectorField& U
)
:
    name_(name),
    rheology_(rheology),
    lambda_
    (
        IOobject
        (
            name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    expLambdaSource_
    (
        IOobject
        (
            "lambdaSource",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("lambdaSource", lambda_.dimensions()/dimTime, 0.0)
    ),
    impLambdaSource_
    (
        fvm::Sp(dimensionedScalar("rtzero", dimless/dimTime, 0.0), lambda_)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::structureModel::correct(const Foam::surfaceScalarField& phi)
{
    bool stabDiff = false;
    dimensionedScalar D = dimensionedScalar("D", dimArea/dimTime, 1.0e-6);

    Info<< "Shear rate min/max = "
        << gMin(rheology_.sr()) << "/"
        << gMax(rheology_.sr()) << endl;
    //Info<< "Explicit lambda source min/max = "
    //    << gMin(expLambdaSource_) << "/"
    //    << gMax(expLambdaSource_) << endl;

    if(lambda_.db().foundObject<solutionControl>("simpleControl"))
    {
        fvScalarMatrix lambdaEqn
        (
            fvm::div(phi, lambda_)
            ==
            expLambdaSource_
            - impLambdaSource_
        );
        if(stabDiff == true)
        {
            lambdaEqn += 
                - fvm::laplacian(D, lambda_, "fvm::laplacian(D,lambda)")
                //- fvc::div(D*dev(T(fvc::grad(lambda_))), "div(D*dev(T(grad(lambda))))")
                + fvc::div(D*fvc::grad(lambda_), "div(D,grad(lambda))");
        }

        lambdaEqn.relax();
        lambdaEqn.solve();
    }
    else
    {
        fvScalarMatrix lambdaEqn
        (
            fvm::ddt(lambda_)
            + fvm::div(phi, lambda_)
            ==
            expLambdaSource_
            - impLambdaSource_
        );

        if(stabDiff == true)
        {
            lambdaEqn += 
                - fvm::laplacian(D, lambda_, "fvm::laplacian(D,lambda)")
                //- fvc::div(D*dev(T(fvc::grad(lambda_))), "div(D*dev(T(grad(lambda))))")
                + fvc::div(D*fvc::grad(lambda_), "div(D,grad(lambda))");
        }

        lambdaEqn.relax();
        lambdaEqn.solve();
    }

    //Info<< "Lambda min/max = "
    //    << gMin(lambda_) << "/"
    //    << gMax(lambda_) << endl;

}

// ************************************************************************* //


/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "myCfdemCloud.H"
#include "explicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"
#include "cellSet.H"
#include "fvIOoptionList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFvOptions.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    myCfdemCloud particleCloud(mesh);
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        particleCloud.clockM().start(1,"Global");

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothen(particleCloud.forceM(0).impParticleForces());
        }
    
        Info << "update f.internalField()" << endl;
        f = particleCloud.momCoupleM(0).expMomSource();
        f.correctBoundaryConditions();

       //Force Checks
       vector fImpTotal(0,0,0);
       vector fExpTotal = sum(mesh.V()*f.internalField());
       reduce(fExpTotal, sumOp<vector>());
       Info << "TotalForceExp: " << fExpTotal << endl;
       Info << "TotalForceImp: " << fImpTotal << endl;

        //#include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
//                + turbulence->divDevReff(U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                    f/rho
                  + fvOptions(U)
                );

                UEqn.relax();
	            fvOptions.constrain(UEqn);

                if (momentumPredictor && (modelType=="B" || modelType=="Bfull"))
		        {
                    solve(UEqn == - fvc::grad(p));
	                fvOptions.correct(U);
		        }
                else if (momentumPredictor)
		        {
                    solve(UEqn == - voidfraction*fvc::grad(p));
	                fvOptions.correct(U);
		        }

                // --- PISO loop

                //for (int corr=0; corr<nCorr; corr++)
                int nCorrSoph = nCorr + 5 * pow((1-particleCloud.dataExchangeM().timeStepFraction()),1);

                for (int corr=0; corr<nCorrSoph; corr++)
                {
                    volScalarField rUA = 1.0/UEqn.A();

                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

                    U = rUA*UEqn.H();

                    #ifdef version23
                    phi = ( fvc::interpolate(U*voidfraction) & mesh.Sf() )
                        + rUAfvoidfraction*fvc::ddtCorr(U, phi);
                    #else
                    phi = ( fvc::interpolate(U*voidfraction) & mesh.Sf() )
                        + fvc::ddtPhiCorr(rUAvoidfraction, U, phi);
                    #endif

                    fvOptions.makeRelative(phi);
                    adjustPhi(phi, U, p);

                    //surfaceScalarField phiS(fvc::interpolate(Us*voidfraction) & mesh.Sf());
                    surfaceScalarField phiGes = phi; //+ rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                    if (modelType=="A")
                        rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                    // Non-orthogonal pressure corrector loop
                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    {
                        // Pressure corrector
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rUAvoidfraction, p) == fvc::div(phiGes) + particleCloud.ddtVoidfraction()
                        );
                        pEqn.setReference(pRefCell, pRefValue);

                        if
                        (
                            corr == nCorr-1
                         && nonOrth == nNonOrthCorr
                        )
                        {
                            pEqn.solve(mesh.solver("pFinal"));
                        }
                        else
                        {
                            pEqn.solve();
                        }

                        if (nonOrth == nNonOrthCorr)
                        {
                            phiGes -= pEqn.flux();
                            phi = phiGes;
                        }

                    } // end non-orthogonal corrector loop

                    #include "continuityErrorPhiPU.H"

                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p); //- Ksl/rho*Us*rUA;
                    else
                        U -= voidfraction*rUA*fvc::grad(p); //- Ksl/rho*Us*rUA;

                    U.correctBoundaryConditions();

                } // end piso loop
            }

            turbulence->correct();
            fvOptions.correct(U);
        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

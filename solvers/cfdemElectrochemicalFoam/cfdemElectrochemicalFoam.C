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

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/
/*Mohsen
#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "rhoChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"*/



#include "fvCFD.H"
//#include "fvm.H"
//#include "fvOption.H"
#include "fvIOoptionList.H"
#include "meshToMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "concChemistryModel.H"
#include "multivariateScheme.H"
//#include "pimpleControl.H"

#include "myCfdemCloud.H"
#include "implicitCouple.H"
#include "smoothingModel.H"
#include "forceModel.H"
// ----------------------------- code addition ----------------------------- //
#include "calculatedGradientFvPatchScalarField.H"
#include "electrochemistryModel.H"
//#include "liquidMultiSpeciesTransportModel.H"
#include "OFstream.H"
#include <iostream>
#include <fstream>
// ------------------------------------------------------------------------- //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// A.Alexiou 2015, to use OF 2.1 solver style define RHO_EQN_2_1
//#define RHO_EQN_2_1


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readChemistryProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFvOptions.H"
    #include "readTimeControls.H"
    //#include "CourantNo.H"
    //#include "setInitialDeltaT.H"

    //pimpleControl pimple(mesh);


    #include "checkModelType.H"

    Info << "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        particleCloud.clockM().start(1,"Global");
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"
        //#include "setDeltaT.H"
        #include "meshToMesh.H"
        //runTime++;

        // do particle stuff
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        //bool hasEvolved = true;

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothen(particleCloud.forceM(0).impParticleForces());
        }
    
        Info << "update Ksl.internalField()" << endl;
        Ksl = particleCloud.momCoupleM(0).impMomSource();
        Ksl.correctBoundaryConditions();

       //Force Checks
       vector fTotal(0,0,0);
       vector fImpTotal = sum(mesh.V()*Ksl.internalField()*(Us.internalField()-U.internalField()));
       reduce(fImpTotal, sumOp<vector>());
       Info << "TotalForceExp: " << fTotal << endl;
       Info << "TotalForceImp: " << fImpTotal << endl;

        #include "solverDebugInfo.H"

        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");

        porosityFactor = pow(voidfraction, BruggemanCoeff);
        // get scalar source from DEM        
        //particleCloud.forceM(1).manipulateScalarField(Tsource);
        //Tsource.correctBoundaryConditions();

        // solve scalar transport equation
        //fvScalarMatrix TEqn
        //(
        //   fvm::ddt(voidfraction,T) - fvm::Sp(fvc::ddt(voidfraction),T)
        // + fvm::div(phi, T) - fvm::Sp(fvc::div(phi),T)
        // - fvm::laplacian(DT*voidfraction, T)
        // ==
        //   Tsource
        //);
        //TEqn.relax();
        //TEqn.solve();
        
        #include "chemistry.H"
        
        //fv::IOoptionList fvOptions(mesh); // for compatibility with OF 2.3, since rhoEqn.H changed

        //#include "rhoEqn.H" // OF header/source

        //while (pimple.loop()) 
        //{

            #include "CEqn.H"
            #include "ElectricEqn.H"

            //volScalarField rhoVoidfraction = rho*voidfraction;
            //#include "hsEqn.H"

            if(particleCloud.solveFlow())
            {
                #include "UEqn.H"
                
                int nCorrSoph = nCorr + 5 * pow((1-particleCloud.dataExchangeM().timeStepFraction()),1);

                for (int corr=0; corr<nCorrSoph; corr++)
                //while (pimple.correct())
                {

                    #include "pEqn.H"
                    
                }
                //if (pimple.turbCorr())
                //{ 
                    turbulence->correct();
                //}
                fvOptions.correct(U);
            }
            else
            {
                Info << "skipping flow solution." << endl;
            }
        //}
        
        scalar iAvg = 0.0;
        scalar phiEsAvg = electrodes[0].potential();
        {
            vectorField ilp = mesh.boundary()["separator"].patchField<volVectorField, vector>(il);
            scalarField Ap = mesh.boundary()["separator"].magSf();
            tmp<vectorField> tn = mesh.boundary()["separator"].nf();
            vectorField n = tn();
            scalarField iNormal = ilp&n;
            label faceCount = 0;
            forAll(iNormal, faceI)
            {
                ++faceCount;
                iAvg += iNormal[faceI];
            }
            iAvg /= faceCount;
            Info<< "Average cell current density = " << iAvg << endl; 
            Info << "Average electrode potential = " << phiEsAvg << endl; 
            if (timeStepCounter > startCounter)
            { 
                iAvgTime += iAvg;
                phiEsAvgTime += phiEsAvg;
            }
            ++timeStepCounter;
        }
        runTime.write();
        if (runTime.write())
        {
            //thermo.T().write();
            //chemistry.dQ()().write();
            //forAll(composition.C(), i)
            //{
            //    composition.C()[i].write();
            //}
            //iTrans.write();
            //is.write();
        }
        
        // Calculate average current density at separator

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    fileName outputFile(runTime.path()+"/CurrentDensity-Potential.txt");
    OFstream os(outputFile);
    os << "Average cell current density: " 
        << iAvgTime/scalar(timeStepCounter-startCounter) << endl;
    os << "Average electrode potential: " 
        << phiEsAvgTime/scalar(timeStepCounter-startCounter) << endl;

    Info << "Average cell current density: " 
        << iAvgTime/scalar(timeStepCounter-startCounter) << endl;
    Info << "Average electrode potential: " 
        << phiEsAvgTime/scalar(timeStepCounter-startCounter) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

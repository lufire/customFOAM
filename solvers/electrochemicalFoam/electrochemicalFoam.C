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
    electrochemicalFoam

Description
    Solver for transport in electrolytes with electrochemical surface reactions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
//#include "turbulentFluidThermoModel.H"
#include "concReactionThermo.H"
#include "pimpleControl.H"
//#include "pressureControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
//#include "turbulenceModel.H"
#include "multivariateScheme.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// ----------------------------- code addition ----------------------------- //
#include "calculatedGradientFvPatchScalarField.H"
#include "electrochemistryModel.H"
#include "OFstream.H"
#include <iostream>
#include <fstream>
// ------------------------------------------------------------------------- //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// A.Alexiou 2015, to use OF 2.1 solver style define RHO_EQN_2_1
//#define RHO_EQN_2_1


int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    //#include "readChemistryProperties.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            //#include "chemistry.H"
            #include "CEqn.H"
            #include "ElectricEqn.H"
            //#include "hsEqn.H"
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        //scalar phiEsAvg = electrodes[0].potential();
        scalar phiEsAvg = 0.0;
        scalar iAvg = 0.0;
        
        {
            vectorField ilp = 
                mesh.boundary()["currentCollector"].patchField<volVectorField, vector>(il);
            tmp<vectorField> tn = mesh.boundary()["currentCollector"].nf();
            const vectorField& n = tn();
            scalarField iNormal = ilp&n;
            const scalarField& magSf = mesh.boundary()["currentCollector"].magSf();
            label faceCount = 0;
            scalar patchArea = 0.0;
            forAll(iNormal, faceI)
            {
                patchArea += magSf[faceI];
                iAvg += iNormal[faceI]*magSf[faceI];

            }
            if(patchArea > 0.0) 
            {
                iAvg /= patchArea;
            }
            //Info << "Average cell current density:  " << iAvg << endl; 
            //Info << "Average electrode potential: " << phiEsAvg << endl; 
            if (timeStepCounter > startCounter)
            { 
                iAvgTime += iAvg;
                //phiEsAvgTime += phiEsAvg;
            }
            ++timeStepCounter;
        }
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    //fileName outputFile(runTime.path()+"/CurrentDensity-Potential.txt");
    //OFstream os(outputFile);
    //os << "Average cell current density: " 
    //    << iAvgTime/scalar(timeStepCounter-startCounter) << endl;
    //os << "Average electrode potential: " 
    //    << phiEsAvgTime/scalar(timeStepCounter-startCounter) << endl;

    //Info << "Average cell current density: " 
    //    << iAvgTime/scalar(timeStepCounter-startCounter) << endl;
    //Info << "Average electrode potential: " 
    //    << phiEsAvgTime/scalar(timeStepCounter-startCounter) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

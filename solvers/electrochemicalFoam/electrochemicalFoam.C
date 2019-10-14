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
#include "fvm.H"
#include "fvOption.H"
#include "fvIOoptionList.H"
#include "meshToMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "concChemistryModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"

// ----------------------------- code addition ----------------------------- //
#include "calculatedGradientFvPatchScalarField.H"
#include "linearGradientFvPatchScalarField.H"
#include "electrochemistryModel.H"
//#include "conductivityModel.H"
//#include "diffusivityModel.H"
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
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);


    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        while (pimple.loop()) 
        {
            //surfaceScalarField phiByRho = phi/fvc::interpolate(rho);

            //tmp<fv::convectionScheme<scalar> > mvConvection
            //(
            //    fv::convectionScheme<scalar>::New
            //    (
            //        mesh,
            //        fields,
            //        phiByRho,
            //        mesh.divScheme("div(phi,Ci_h)")
            //    )
            //);

            #include "chemistry.H"
            label iter = 1;
            label iterMax = 10; 
            scalar alpha = 0.0;
            while(iter <= iterMax)
            {
                //if(iter == iterMax)
                //{
                //    alpha = 1.0;
                //}
                alpha = Foam::pow(scalar(iter)/scalar(iterMax), 2.0);
                #include "CEqn.H"
                #include "ElectricEqn.H"
                iter++;
            }
            
            //#include "hsEqn.H"
            #include "UEqn.H"
            while (pimple.correct())
            {

                #include "pEqn.H"
                
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        //scalar phiEsAvg = electrodes[0].potential();
        scalar phiEsAvg = 0.0;
        scalar iAvg = 0.0;
        
        {
            vectorField ilp = mesh.boundary()["currentCollector"].patchField<volVectorField, vector>(il);
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
            Info << "Average cell current density = " << iAvg << endl; 
            //Info << "Average electrode potential = " << phiEsAvg << endl; 
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

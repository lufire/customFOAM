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
//#include "meshToMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "concChemistryModel.H"
#include "fvIOoptionList.H"
#include "pimpleControl.H"
//#include "fixedFluxPressureFvPatchScalarField.H"

// ----------------------------- code addition ----------------------------- //
#include "calculatedGradientFvPatchScalarField.H"
#include "electrochemistryModel.H"
// ------------------------------------------------------------------------- //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


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
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);


    Info << "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop()) 
        {

            #include "chemistry.H"
            //#include "hsEqn.H"
            #include "UEqn.H"
            #include "CEqn.H"
            #include "ElectricEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {

                #include "pEqn.H"
                
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

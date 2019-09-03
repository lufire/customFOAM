/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::boundaryPatchElectrode

Description
    Base-class for all electrochemistry models.

SourceFiles
    boundaryPatchElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryPatchElectrode_H
#define boundaryPatchElectrode_H

#include "electrodeModel.H"
#include "resistanceModel.H"
//#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

/*---------------------------------------------------------------------------*\
                Class boundaryPatchElectrode Declaration
\*---------------------------------------------------------------------------*/

class boundaryPatchElectrode
:
    public electrodeModel
{
    // Private data
  
        //- Electrode potential
        scalarField phiEs_;
        
        //- Electrode over potential
        scalarField eta_;
        
        //- Electrode potential
        scalar avgPhiEs_;
        
        //- Target current density 
        scalar currD_;
        
        //- Resistance model
        autoPtr<resistanceModel> resistanceModel_;

        bool fluxCorrector_;

    // Private Member Functions

        //- Disallow copy construct
        boundaryPatchElectrode(const boundaryPatchElectrode&);

        //- Disallow default bitwise assignment
        void operator=(const boundaryPatchElectrode&);

        //- Initialize members variables
        void initializeMembers();

        //- Correct the electrode surface area per cell volume
        void correctElectrodeSurface();

        //- Correct the electrode surface area per cell volume
        void correctElectrodeConductivity();

        //- Correct the transfer current density
        void correctTransferCurrentDensity();

        //- Correct the transfer current density
        void correctPotential();

public:

    //- Runtime type information
    TypeName("boundaryPatchElectrode");


    // Constructors

        //- Construct from components
        boundaryPatchElectrode
        (
            const dictionary& dict, 
            electrolyteModels::electrolyteModel& electrolyte 
        );


    //- Destructor
    virtual ~boundaryPatchElectrode();


    // Member Functions

        //- Update electrochemistry model properties
        void correct(){};

        virtual void correctSpeciesFlux();

        virtual void correctElectricity();

        virtual void correctSpecies(PtrList<volScalarField>& C);


        inline scalar potential()
        {
            return avgPhiEs_;
        };
        //- Read transportProperties dictionary
        //bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "boundaryPatchElectrodeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

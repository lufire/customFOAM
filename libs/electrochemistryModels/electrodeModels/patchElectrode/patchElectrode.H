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
    Foam::patchElectrode

Description
    Base-class for all electrochemistry models.

SourceFiles
    patchElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef patchElectrode_H
#define patchElectrode_H

#include "electrodeModel.H"
//#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

/*---------------------------------------------------------------------------*\
                Class patchElectrode Declaration
\*---------------------------------------------------------------------------*/

class patchElectrode
:
    public electrodeModel
{
    // Private data
  
        //- Electrode potential
        dimensionedScalar phiEs_;

    // Private Member Functions

        //- Disallow copy construct
        patchElectrode(const patchElectrode&);

        //- Disallow default bitwise assignment
        void operator=(const patchElectrode&);

        //- Correct the transfer current density
        void correctTransferCurrentDensity();

        //- Correct the electrode surface area per cell volume
        void correctElectrodeSurface();

        //- Correct the electrode surface area per cell volume
        void correctElectrodeConductivity();


public:

    //- Runtime type information
    TypeName("patchElectrode");


    // Constructors

        //- Construct from components
        patchElectrode
        (
            const dictionary& dict, 
            electrolyteModels::electrolyteModel& electrolyte
        );


    //- Destructor
    virtual ~patchElectrode();


    // Member Functions

        //- Return average electrode potential
        virtual scalar potential()
        {
            return 0.0;
        };

        //- Update electrochemistry model properties
        void correct();

        //- Read transportProperties dictionary
        //bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "patchElectrodeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

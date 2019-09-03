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
    Foam::cfdemPorousElectrode

Description
    Base-class for all electrochemistry models.

SourceFiles
    cfdemPorousElectrode.C

\*---------------------------------------------------------------------------*/

#ifndef cfdemPorousElectrode_H
#define cfdemPorousElectrode_H

#include "electrodeModel.H"
//#include "myCfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{

/*---------------------------------------------------------------------------*\
                Class cfdemPorousElectrode Declaration
\*---------------------------------------------------------------------------*/

class cfdemPorousElectrode
:
    public electrodeModel
{
    // Private data
  
        //- CFDEM particle cloud data
        const myCfdemCloud& particleCloud_;

        //- Percolation field 
        volScalarField percolationField_;

        //- Electric potential field of electrode
        volScalarField phiEs_;

        //- Interval for calculation of the percolation field
        const label percolationInterval_;

        //- Counter for tracking the interval
        label intervalCounter_;

        //- Counter for tracking the interval
        scalar contactDistance_;

    // Private Member Functions

        //- Disallow copy construct
        cfdemPorousElectrode(const cfdemPorousElectrode&);

        //- Disallow default bitwise assignment
        void operator=(const cfdemPorousElectrode&);

        //- Correct the transfer current density
        void correctTransferCurrentDensity();

        //- Correct the electrode surface area per cell volume
        void correctPercolation();

        //- Correct the electrode surface area per cell volume
        void correctElectrodeSurface();

        //- Correct the electrode surface area per cell volume
        void correctElectrodeConductivity();


public:

    //- Runtime type information
    TypeName("cfdemPorousElectrode");


    // Constructors

        //- Construct from components
        cfdemPorousElectrode
        (
            const dictionary& dict, 
            electrolyteModels::electrolyteModel& electrolyte,
            const myCfdemCloud& particleCloud
        );


    //- Destructor
    virtual ~cfdemPorousElectrode();


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

//#include "cfdemPorousElectrodeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

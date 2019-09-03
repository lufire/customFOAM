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
    Foam::electrochemistryModel

Description
    Base-class for all electrochemistry models.

SourceFiles
    electrochemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef electrochemistryModel_H
#define electrochemistryModel_H

#include "primitiveFieldsFwd.H"
#include "volFields.H"
//#include "surfaceFieldsFwd.H"
#include "fvMesh.H"
#include "IOdictionary.H"
//#include "typeInfo.H"
//#include "myCfdemCloud.H"
#include "concReactionThermo.H"
#include "hashedWordList.H"
#include "electrodeModel.H"
#include "electrolyteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class myCfdemCloud;

namespace electrochemistryModels 
{



/*---------------------------------------------------------------------------*\
                Class electrochemistryModel Declaration
\*---------------------------------------------------------------------------*/

class electrochemistryModel
{
    // Private data
  
        //- Mesh
        const fvMesh& mesh_;

        //- Electrochemistry dictionary
        const IOdictionary electrochemDict_;

        //- Bool to set chemistry model active 
        bool isActive_;

        //- Bool to set chemistry model active 
        bool currentControl_;

        //- Electrolyte model
        autoPtr<electrolyteModels::electrolyteModel> electrolyte_;

        //- List of electrode names
        const hashedWordList electrodeTable_;
    
        //- Electrode model pointer list
        PtrList<electrodeModels::electrodeModel> electrodes_;

    // Private Member Functions

        //- Disallow copy construct
        electrochemistryModel(const electrochemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const electrochemistryModel&);


public:

    //- Runtime type information
    TypeName("electrochemistryModel");

    // Constructors

        //- Construct from components
        electrochemistryModel
        (
            const concReactionThermo& thermo,
            const myCfdemCloud& particleCloud
        );

        electrochemistryModel
        (
            const concReactionThermo& thermo
        );


    //- Destructor
    virtual ~electrochemistryModel();


    // Member Functions

        //- Return the const mass-fraction field for a specie given by index
        inline const electrolyteModels::electrolyteModel& electrolyte() const;
        
        //- Return the const mass-fraction field for a specie given by index
        inline electrolyteModels::electrolyteModel& electrolyte();
        
        //- Return the mass-fraction fields
        inline PtrList<electrodeModels::electrodeModel>& electrodes();

        //- Return the const mass-fraction fields
        inline const PtrList<electrodeModels::electrodeModel>& electrodes() const;

        //- Return the mass-fraction field for a specie given by index
        inline electrodeModels::electrodeModel& electrode(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const electrodeModels::electrodeModel& electrode(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline electrodeModels::electrodeModel& electrode(const word& electrodeName);

        //- Return the const mass-fraction field for a specie given by name
        inline const electrodeModels::electrodeModel& electrode(const word& electrodeName) const;

        //- Return the const mass-fraction field for a specie given by name
        inline bool active();

        //- Update electrochemistry model properties
        void correct();

        //- Read transportProperties dictionary
        //bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "electrochemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

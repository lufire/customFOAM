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

Class
    Foam::resistanceModel

Description
    Abstract base class for diffusivity calculation

SourceFiles
    resistanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef resistanceModel_H
#define resistanceModel_H

//#include "specie.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
//#include "scalarField.H"
//#include "volFieldsFwd.H"
//#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
//#include "tmp.H"
#include "autoPtr.H"
//#include "molecularWeights.H"
//#include "porosityModelList.H"
#include "volFields.H"
#include "concReactionThermo.H"
//#include "electrodeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrodeModels
{
//- Forward declaration of classes
class electrodeModel;

/*---------------------------------------------------------------------------*\
                    Class resistanceModel Declaration
\*---------------------------------------------------------------------------*/

class resistanceModel
{

protected:

    // Protected data
    
        //- Transport propertied dictionary
        const dictionary dict_;

        //- Electrode object 
        const electrodeModel& electrode_;

        //- Thermo object
        const concReactionThermo& thermo_;

        //- Mesh
        const fvMesh& mesh_;

        //- Species table
        const speciesTable& species_;

        //- Concentration at electrode surface
        PtrList<scalarField> Cs_;

        //- Electric potential of electrolyte at surface
        scalarField phiEls_;

        //- Rotation matrix to convert vectors to local coordinate system
        //- based on patch normal vector (i.e. local x-coordinate)
        //tensorField rotationMatrix_;


    // Private Member Functions

        //- Disallow copy construct
        resistanceModel(const resistanceModel&);

        //- Disallow default bitwise assignment
        void operator=(const resistanceModel&);


public:

    //- Runtime type information
    TypeName("resistanceModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            resistanceModel,
            dictionary,
            (
                const dictionary& dict,
                const electrodeModel& electrode
            ),
            (dict, electrode)
        );

    // Selectors

        static autoPtr<resistanceModel> New
        (
            const dictionary& dict,
            const electrodeModel& electrode 
        );

    // Constructors

        //- Construct from components
        resistanceModel
        (
            const dictionary& dict,
            const electrodeModel& electrode 
        );

    // Destructor

        virtual ~resistanceModel()
        {}

    // Member Functions
    //
        //- Correct species concentrations 
        virtual void correctSpecies() = 0;

        //- Correct electric potential 
        virtual void correctPotential() = 0;

        //- Correct the species flux 
        virtual void correctSpeciesFlux(vectorField& N) = 0;

        //- Calculate diffusivity 
        //virtual void update() = 0;

        //- Return const species concentrations at electrode surface
        inline const PtrList<scalarField>& Cs() const;

        //- Return const species concentrations at electrode surface
        inline PtrList<scalarField>& Cs();

        //- Return const species concentrations at electrode surface
        inline const scalarField& Cs(const label i) const;

        //- Return const species concentrations at electrode surface
        inline scalarField& Cs(const label i);

        //- Return const electric potential of electrolyte at electrode surface
        inline const scalarField& phiEls() const;

        //- Return electric potential of electrolyte at electrode surface
        inline scalarField& phiEls();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrodeModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
#include "resistanceModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

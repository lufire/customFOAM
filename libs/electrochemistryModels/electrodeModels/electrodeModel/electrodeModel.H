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
    Foam::electrodeModel

Description
    Abstract base-class for all electrode models.

SourceFiles
    electrodeModel.C

\*---------------------------------------------------------------------------*/

#ifndef electrodeModel_H
#define electrodeModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
//#include "surfaceFieldsFwd.H"
#include "fvMesh.H"
#include "dictionary.H"
//#include "typeInfo.H"
#include "concReactionThermo.H"
#include "runTimeSelectionTables.H"
#include "electrolyteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class myCfdemCloud;

namespace electrochemistryModels 
{
namespace electrodeModels
{

/*---------------------------------------------------------------------------*\
                Class electrodeModel Declaration
\*---------------------------------------------------------------------------*/

class electrodeModel
{

protected:

    // Protected data
  
        //- Electrode subdictionary
        const dictionary& electrodeDict_;

        //- Reference to electrolyte object 
        electrolyteModels::electrolyteModel& electrolyte_;

        //- Mesh
        const fvMesh& mesh_;

        //- Patch name 
        const word patchName_;

        //- Patch ID
        const label patchID_;

        //- Bool to set control mode 
        bool currentControl_;

        //- Electrode surface density (volume specific)
        //volScalarField surfaceDensity_;
        scalarField surfaceDensity_;

        //- Temperature field 
        const volScalarField& T_;

        //- Oxidant species field name
        const word oxName_;
        
        //- Reductant species field name
        const word redName_;
        
        //- Species concentrations 
        const PtrList<volScalarField>& C_;

        //- Number of transferred electrons
        const label electronNumber_;
 
        //- Stoichiometric coefficients for electrochemical reaction
        PtrList<label> stoichCoeff_;

        //- Reaction orders for electrochemical reaction
        PtrList<label> reactionOrder_;

        //- Species refence concentrations 
        PtrList<scalar> CRef_;
        
        //- Exchange current density
        const scalar iEx_;
        
        //- Equilibrium potential (i.e. open circuit voltage)
        const scalar eqPotential_;
        
        //- Anodic charge transfer coefficient
        const scalar alphaA_;
        
        //- Cathodic charge transfer coefficient
        const scalar alphaC_;

        //- Charge transfer current density (volumetric)
        volScalarField iTrans_;

        //- Current density field of electrode
        volVectorField is_;
         
        //- Index of oxidant species
        label oxID_;
        
        //- Index of reductant species
        label redID_;

        //- Electrode conductivity
        const dimensionedScalar sigma_;

        //- Effective electrode conductivity
        volScalarField sigmaEff_;

private:

    // Private Member Functions

        //- Disallow copy construct
        electrodeModel(const electrodeModel&);

        //- Disallow default bitwise assignment
        void operator=(const electrodeModel&);

        void initializeMembers();

protected:

    // Protected member functions

        //- Correct the transfer current density
        virtual void correctTransferCurrentDensity() = 0;

        //- Correct the electrode surface area per cell volume
        virtual void correctElectrodeSurface() = 0;

        //- Correct the electrode surface area per cell volume
        virtual void correctElectrodeConductivity() = 0;


public:

    //- Runtime type information
    TypeName("electrodeModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            electrodeModel,
            electrolyte,
            (
                const dictionary& dict,
                electrolyteModels::electrolyteModel& electrolyte 
            ),
            (dict, electrolyte)
        );

        //declareRunTimeSelectionTable
        //(
        //    autoPtr,
        //    electrodeModel,
        //    dictionary,
        //    (
        //        const dictionary& electrodeDict,
        //        const concReactionThermo& thermo
        //    ),
        //    (electrodeDict, thermo)
        //);

        declareRunTimeSelectionTable
        (
            autoPtr,
            electrodeModel,
            myCfdemCloud,
            (
                const dictionary& dict,
                electrolyteModels::electrolyteModel& electrolyte,
                const myCfdemCloud& particleCloud
            ),
            (dict, electrolyte, particleCloud)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<electrodeModel> New
        (
            const dictionary& dict,
            electrolyteModels::electrolyteModel& electrolyte 
        );

        //static autoPtr<electrodeModel> New
        //(
        //    const dictionary& electrodeDict,
        //    const concReactionThermo& thermo
        //);

        static autoPtr<electrodeModel> New
        (
            const dictionary& dict,
            electrolyteModels::electrolyteModel& electrolyte,
            const myCfdemCloud& particleCloud
        );

    // Constructors

        //- Construct from components
        electrodeModel
        (
            const dictionary& dict, 
            electrolyteModels::electrolyteModel& electrolyte
        );

        //electrodeModel
        //(
        //    const dictionary& electrodeDict, 
        //    const concReactionThermo& thermo
        //);

        electrodeModel
        (
            const dictionary& dict, 
            electrolyteModels::electrolyteModel& electrolyte,
            const myCfdemCloud& particleCloud
        );


    //- Destructor
    virtual ~electrodeModel();


    // Member Functions

        //- Return electrolyte model 
        inline const electrolyteModels::electrolyteModel& electrolyte() const; 

        //- Return transfer current density
        inline volScalarField& iTrans(); 

        //- Return electric potential field of electrode
        //inline volScalarField& phiEs();

        //- Return electronic current density
        inline volVectorField& is();

        //- Return stoichiometric coefficients
        inline const PtrList<label>& stoichCoeff() const;

        //- Return stoichiometric coefficient of specie
        inline label stoichCoeff(const label i) const;

        //- Return stoichiometric coefficient of specie
        inline label electronNumber() const;

        //- Return stoichiometric coefficient of specie
        inline label oxidantIndex() const;

        //- Return stoichiometric coefficient of specie
        inline label reductantIndex() const;

        //- Return stoichiometric coefficient of specie
        inline label patchIndex() const;

        //- Return stoichiometric coefficient of specie
        inline const word& patchName() const;

        //- Return stoichiometric coefficient of specie
        inline const volScalarField& sigmaEff() const;

        virtual void correctSpeciesFlux(){};

        virtual void correctElectricity(){};

        //virtual void correctSpecies(){};
        virtual void correctSpecies(PtrList<volScalarField>& C){};

        //- Return average electrode potential
        virtual scalar potential() = 0;

        //- Update electrochemistry model properties
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        //bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrodeModels
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "electrodeModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

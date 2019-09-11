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
    Foam::liquidMultiSpeciesTransportModel

Description
    Abstract base class for liquidMultiSpecies mass transport model models.

SourceFiles
    liquidMultiSpeciesTransportModelI.H
    liquidMultiSpeciesTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef liquidMultiSpeciesTransportModel_H
#define liquidMultiSpeciesTransportModel_H

#include "fvc.H"

#include "concChemistryModel.H"
#include "turbulenceModel.H"
#include "diffusivityModel.H"
#include "concMultiComponentMixture.H"
//#include "concReactingMixture.H"
#include "porosityModelList.H"
#include "concReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class liquidMultiSpeciesTransportModel Declaration
\*---------------------------------------------------------------------------*/

class liquidMultiSpeciesTransportModel
:

    public IOdictionary
{
    
protected:  
  
    // Protected data

        //- Mesh
        const fvMesh& mesh_;

        // Porous zones
        porosityModelList pZones_ ;

        // Reference to the thermodynamic package
        concReactionThermo& thermo_;

        //- Index of onert species
        const scalar inertIndex_;

        //- Reference to the turbulence model
        const incompressible::turbulenceModel& turbulence_;

        //- Molar fractions [-]
        PtrList<volScalarField> x_;

        //- Species mass fluxes [kg/s]
        PtrList<surfaceScalarField> n_;

        //- Species rates [kg/m^3/s]
        PtrList<volScalarField> Sy_;

        //- Diffusivity models 
        autoPtr<diffusivityModel> DijModel_;

        //- liquidMultiSpecies diffusion coefficients
        PtrList<volScalarField> D_;
        
        //- Valence of the ionic species
        PtrList<label> z_;            


private:
  
    // Private member functions
      
        //- Disallow copy construct
        liquidMultiSpeciesTransportModel(const liquidMultiSpeciesTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const liquidMultiSpeciesTransportModel&);

//        Foam::liquidMultiSpeciesTransportModel::RR
//        (
//            const word Yname,
//            const concChemistryModel& chemistry,
//            const label i
//        ) const


public:
 
    //- Runtime type information
    TypeName("liquidMultiSpeciesTransportModel");

      
    // Declare run-time constructor selection table
       
#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        liquidMultiSpeciesTransportModel,
        fvMesh,
        (
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        ),
        (thermo, turbulence)
    );
#endif
  

    // Constructors	  

        //- Construct from components
        liquidMultiSpeciesTransportModel
        (
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        );


    // Selectors	  
  
        //- Return a reference to the selected liquidMultiSpecies transport model

        static autoPtr<liquidMultiSpeciesTransportModel> New
        (
            concReactionThermo& thermo,
            const incompressible::turbulenceModel& turbulence
        );


    //- Destructor

        virtual ~liquidMultiSpeciesTransportModel()
        {}


    // Member Functions

        // Access   
           
            //- Return the const binary diffusion coefficient field
            inline const volScalarField& Dij(label i, label j) const;
            
            //- Return the const diffusion coefficient field
            inline const volScalarField& Di(label i) const;
            
            //- Return the const list of diffusion coefficient fields
            inline const PtrList<volScalarField>& D() const;
            
            //- Return the porous zones
            inline porosityModelList& pZones();
            
            //- Return the table of species
            inline const speciesTable& species() const;

            //- Return the const total mass flux field for a species
            inline const surfaceScalarField& n(const label i) const;
            
            //- Return the const valence of ionic species
            inline const label& z(label i) const;                       

            //- Return the const list of diffusion coefficient fields
            inline const PtrList<label>& z() const;

            //- Return the diffusive mass flux field for a specie given by index
            surfaceScalarField j(const label i);

            //- Molar weight for a specie given by index
            virtual scalar W(label i) const = 0;

            //- Specific enthalpy for a specie given by index and temperature
            virtual scalar hs(label i, label T) const = 0;


        //- Source term in enthalpy equation due to diffusive fluxes
        //tmp<volScalarField> multiSpeciesHeatSource();

        //- Costruct molar fraction fields
        void updateMolarFractions();         

        //- Update the coefficients of the species equations
        virtual void updateCoefficients() = 0;                

        //- Solve the mass transport equantions without reactions
        virtual scalar correct
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        ) = 0;   

        //- Solve the mass transport equantions with reactions
        virtual scalar correct
        (
            const PtrList<volScalarField>& Y,
            const volScalarField& kappa,
            const concChemistryModel& chemistry,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );
        
        /*
        tmp<Foam::volScalarField> RR
        (
            const word& Yname,
            const concChemistryModel& chemistry,
            const label& i
        );
        */         
        
        ///*
        tmp<Foam::volScalarField> RR
        (
            const word& Yname,
            const concChemistryModel& chemistry,
            const label& i
        ) const;
        //*/    
                
        //- Read transportProperties dictionary
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liquidMultiSpeciesTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

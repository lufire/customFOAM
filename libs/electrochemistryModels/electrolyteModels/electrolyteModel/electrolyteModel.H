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
    Foam::electrolyteModel

Description
    Abstract base class for electrolyte transport model models.

SourceFiles
    electrolyteModelI.H
    electrolyteModel.C

\*---------------------------------------------------------------------------*/

#ifndef electrolyteModel_H
#define electrolyteModel_H

#include "fvc.H"

//#include "concChemistryModel.H"
//#include "turbulenceModel.H"
//#include "diffusivityModel.H"
//#include "conductivityModel.H"
#include "concMultiComponentMixture.H"
//#include "concReactingMixture.H"
#include "porosityModelList.H"
#include "concReactionThermo.H"
#include "diffusivityModel.H"
#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels 
{

// Forward declaration of classes
//class diffusivityModel;
//class conductivityModel;
  
/*---------------------------------------------------------------------------*\
                    Class electrolyteModel Declaration
\*---------------------------------------------------------------------------*/

class electrolyteModel
{
    
protected:  
  
    // Protected data

        //- Electrochemistry dictionary
        const dictionary& electrolyteDict_;

        //- Mesh
        const fvMesh& mesh_;

        // Porous zones
        porosityModelList pZones_ ;

        // Reference to the thermodynamic package
        const concReactionThermo& thermo_;

        //- Species table
        const speciesTable& species_;

        //- Index of inert species
        const word balanceSpecies_;

        //- Index of inert species
        const word inertSpecies_;

        //- Standard concentration 
        const dimensionedScalar Cstd_;

        //- Relative permittivity 
        const dimensionedScalar epsilon_;

        //- Molar fractions [-]
        PtrList<volScalarField> x_;

        //- Species molar fluxes [mol/m^2-s]
        PtrList<volVectorField> N_;

        //- Species rates [mol/m^3/s]
        PtrList<volScalarField> S_;

        //- Valence of the ionic species
        PtrList<label> z_;            

        //- Diffusivity models 
        autoPtr<diffusivityModel> diffModel_;

        //- Ionic conductivity model 
        autoPtr<conductivityModel> kappaModel_;

        //- Electric potential field of electrolyte 
        volScalarField phiE_;

        //- Electric current density [A/m^2]
        volVectorField i_;

        //- Contribution of concentration gradients to 
        //  electric current density [A/m^2]
        volVectorField iD_;

        //- Reference to diffusion coefficients
        //PtrList<volScalarField>& D_;

        //- Pointer to ionic conductivity 
        //tmp<volScalarField> kappa_;


    // Protected member functions
    
        //- Update molar fraction fields
        void updateMolarFractions();         

        //- Update the coefficients of the species equations
        virtual void updateCoefficients() = 0;

private:
  
    // Private member functions
      
        //- Disallow copy construct
        electrolyteModel(const electrolyteModel&);

        //- Disallow default bitwise assignment
        void operator=(const electrolyteModel&);

//        Foam::electrolyteModel::RR
//        (
//            const word Yname,
//            const concChemistryModel& chemistry,
//            const label i
//        ) const


public:
 
    //- Runtime type information
    TypeName("electrolyteModel");

      
    // Declare run-time constructor selection table
       
    declareRunTimeSelectionTable
    (
        autoPtr,
        electrolyteModel,
        dictionary,
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        ),
        (dict, thermo)
    );
  

    // Constructors	  

        //- Construct from components
        electrolyteModel
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );


    // Selectors	  
  
        //- Return a reference to the selected liquidMultiSpecies transport model

        static autoPtr<electrolyteModel> New
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );


    //- Destructor

        virtual ~electrolyteModel()
        {}


    // Member Functions

        // Access   
           
            //- Return the const diffusion coefficient field
            inline const volScalarField& D(label i) const;
            
            //- Return the const list of diffusion coefficient fields
            inline const PtrList<volScalarField>& D() const;

            //- Return the porous zones
            inline porosityModelList& pZones();

            //- Return the thermo model 
            inline const concReactionThermo& thermo() const;
            
            //- Return the table of species
            inline const speciesTable& species() const;

            //- Return the const molar flux field for a species
            inline const PtrList<volVectorField>& N() const;
            
            //- Return the const molar flux field for a species
            inline PtrList<volVectorField>& N();
            
            //- Return the const molar flux field for a species
            inline const volVectorField& N(const label i) const;
            
            //- Return the const molar flux field for a species
            inline volVectorField& N(const label i);
            
            //- Return the const charge number of a species 
            inline const label& z(label i) const;                       

            //- Return the const list of charge numbers
            inline const PtrList<label>& z() const;

            //- Return index for inert species 
            inline const word& balanceSpecies() const;

            //- Return index for inert species 
            inline const word& inertSpecies() const;

            //- Return the const conductivity field 
            inline const volScalarField& kappa() const;
            
            //- Return the standard concentration 
            inline const dimensionedScalar& Cstd() const;
            
            //- Return the relative permittivity 
            inline const dimensionedScalar& epsilon() const;
            
            //- Return reference to the diffusivity model
            inline const diffusivityModel& diffModel() const;
            
            //- Return reference to the conductivity model
            inline const conductivityModel& kappaModel() const;

            //- Return const reference to electric potential field
            inline const volScalarField& phiE() const;

            //- Return reference to electric potential field
            inline volScalarField& phiE();

            //- Return const reference to current density field
            inline const volVectorField& i() const;

            //- Return reference to current density field 
            inline volVectorField& i();

            //- Return const reference to diffusion current density field 
            inline const volVectorField& iD() const;

            //- Return reference to diffusion current density field 
            inline volVectorField& iD();


            ////- Return the diffusive mass flux field for a species given by index
            //surfaceScalarField j(const label i);

            //- Molar weight for a species given by index
            virtual scalar W(label i) const = 0;

            //- Specific enthalpy for a species given by index and temperature
            virtual scalar hs(label i, label T) const = 0;


            //- Source term in enthalpy equation due to diffusive fluxes
            //tmp<volScalarField> multiSpeciesHeatSource();


            ////- Update the coefficients of the species equations
            //virtual void updateCoefficients() = 0;                

            //- Update electrolyte properties 
            virtual void correct() = 0;   

            //virtual scalar correct
            //(
            //    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
            //) = 0;   

        //- Solve the mass transport equantions with reactions
        //virtual scalar correct
        //(
        //    const PtrList<volScalarField>& Y,
        //    const volScalarField& gamma,
        //    const concChemistryModel& chemistry,
        //    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        //);
        
        /*
        tmp<Foam::volScalarField> RR
        (
            const word& Yname,
            const concChemistryModel& chemistry,
            const label& i
        );
        */         
        
        ///*
        //tmp<Foam::volScalarField> RR
        //(
        //    const word& Yname,
        //    const concChemistryModel& chemistry,
        //    const label& i
        //) const;
        //*/    
                
        //- Read transportProperties dictionary
        //virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "electrolyteModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::transferenceNumberModel

Description
    Conductivity model for multispecies binary diffusion coefficients 

SourceFiles
    transferenceNumberModel.C

\*---------------------------------------------------------------------------*/

#ifndef transferenceNumberModel_H
#define transferenceNumberModel_H

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
//#include "electrolyteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

// Forward declaration of classes
class electrolyteModel;

/*---------------------------------------------------------------------------*\
                    Class transferenceNumberModel Declaration
\*---------------------------------------------------------------------------*/

class transferenceNumberModel
{

protected:

    // Protected data
    
        //- Transport propertied dictionary
        const dictionary dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Species table
        const speciesTable& species_;

        //- List of transference numbers 
        PtrList<volScalarField> t_; 

    // Private Member Functions

        //- Disallow copy construct
        transferenceNumberModel(const transferenceNumberModel&);

        //- Disallow default bitwise assignment
        void operator=(const transferenceNumberModel&);


public:

    //- Runtime type information
    TypeName("transferenceNumberModel");

    // Declare run-time constructor selection table
    //
        declareRunTimeSelectionTable
        (
            autoPtr,
            transferenceNumberModel,
            dictionary,
            (
                const dictionary& dict,
                const electrolyteModel& electrolyte
            ),
            (dict, electrolyte)
        )
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            transferenceNumberModel,
            constant,
            (
                const dictionary& dict,
                const concReactionThermo& thermo
            ),
            (dict, thermo)
        );

        //declareRunTimeSelectionTable
        //(
        //    autoPtr,
        //    transferenceNumberModel,
        //    dictionary,
        //    (
        //        const dictionary& dict,
        //        const concReactionThermo& thermo,
        //        const PtrList<label>& z
        //    ),
        //    (dict, thermo, z)
        //);

    // Selectors

        static autoPtr<transferenceNumberModel> New
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte 
        );

        static autoPtr<transferenceNumberModel> New
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );

        //static autoPtr<transferenceNumberModel> New
        //(
        //    const dictionary& dict,
        //    const concReactionThermo& thermo,
        //    const PtrList<label>& z
        //);

    // Constructors

        //- Construct from components
        transferenceNumberModel
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte 
        );

        transferenceNumberModel
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );

        //transferenceNumberModel
        //(
        //    const dictionary& dict,
        //    const concReactionThermo& thermo,
        //    const PtrList<label>& z
        //);

    // Destructor

        virtual ~transferenceNumberModel()
        {}

    // Member Functions

        // Return transference numbers 
        inline const PtrList<volScalarField>&  t() const
        {
            return t_;
        }
        //virtual const PtrList<volScalarField>& t() const;

        //- Calculate conductivity 
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

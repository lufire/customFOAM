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
    Foam::diffusivityModel

Description
    Abstract base class for diffusivity calculation

SourceFiles
    diffusivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusivityModel_H
#define diffusivityModel_H

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
//- Forward declaration of classes
class electrolyteModel;

/*---------------------------------------------------------------------------*\
                    Class diffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class diffusivityModel
{

protected:

    // Protected data
    
        //- Transport propertied dictionary
        const dictionary dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Thermo object
        const concReactionThermo& thermo_;

        //- Species table
        const speciesTable& species_;

        //- Species diffusivity
        PtrList<volScalarField> D_;


    // Private Member Functions

        //- Disallow copy construct
        diffusivityModel(const diffusivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffusivityModel&);


public:

    //- Runtime type information
    TypeName("diffusivityModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            diffusivityModel,
            dictionary,
            (
                const dictionary& dict,
                const electrolyteModel& electrolyte
            ),
            (dict, electrolyte)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            diffusivityModel,
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
        //    diffusivityModel,
        //    dictionary,
        //    (
        //        const dictionary& dict,
        //        const concReactionThermo& thermo,
        //        const dimensionedScalar& Cstd 
        //    ),
        //    (dict, thermo, Cstd)
        //);

    // Selectors

        static autoPtr<diffusivityModel> New
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte 
        );

        static autoPtr<diffusivityModel> New
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );

        //static autoPtr<diffusivityModel> New
        //(
        //    const dictionary& dict,
        //    const concReactionThermo& thermo,
        //    const dimensionedScalar& Cstd 
        //);

    // Constructors

        //- Construct from components
        diffusivityModel
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte 
        );

        diffusivityModel
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );

        //diffusivityModel
        //(
        //    const dictionary& dict,
        //    const concReactionThermo& thermo,
        //    const dimensionedScalar& Cstd 
        //);

    // Destructor

        virtual ~diffusivityModel()
        {}

    // Member Functions

        inline const PtrList<volScalarField>& D() const
        {
            return D_;
        }

        inline const volScalarField& D(label i) const
        {
            return D_[i];
        }
        //- Calculate diffusivity 
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

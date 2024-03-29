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

Namespace
    Foam::binaryDiffusivityModels

Description
    Namespace for binary diffusivity models.

Class
    Foam::binaryDiffusivityModel

Description
    An abstract base class for binary diffusivity models.

SourceFiles
    binaryDiffusivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusivityModel_H
#define binaryDiffusivityModel_H

#include "specie.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "molecularWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class binaryDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class binaryDiffusivityModel
{

protected:

    // Protected data
    
        //- Transport properties dictionary
        const dictionary& dic_;

        //- Species name
        const word& name1_, name2_;

        // Pressure [Pa]
        const volScalarField& p_;

        // Temperature [K]
        const volScalarField& T_;

        // Valence of ionic species [-]
        //label z1_, z2_;
	    


    // Private Member Functions

        //- Disallow copy construct
        binaryDiffusivityModel(const binaryDiffusivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const binaryDiffusivityModel&);


public:

    //- Runtime type information
    TypeName("binaryDiffusivityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binaryDiffusivityModel,
            dictionary,
            (
                const word& name1,
                const word& name2,
                const dictionary& dic,
                const volScalarField& p,
                const volScalarField& T
            ),
            (name1, name2, dic, p, T)
        );


    // Selectors

        //- Return a reference to the selected diffusivity model
        static autoPtr<binaryDiffusivityModel> New
        (
            const word& name1,
            const word& name2,
            const dictionary& dic,
            const volScalarField& p,
           const volScalarField& T
        );


    // Constructors

        //- Construct from components
        binaryDiffusivityModel
        (
            const word& name1,
            const word& name2,
            const dictionary& dic,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        virtual ~binaryDiffusivityModel()
        {}


    // Member Functions
    
        //- Binary diffusion coefficient for patch [m^2/s]
        virtual tmp<scalarField> D
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

         //- Binary diffusion coefficient [m^2/s]
         virtual tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

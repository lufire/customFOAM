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
    Foam::binaryConductivityModels::varyingDiff

Description
    Constant binary diffusivity model.

SourceFiles
    varyingDiff.C

\*---------------------------------------------------------------------------*/

#ifndef varyingDiff_H
#define varyingDiff_H

#include "diffusivityModel.H"
#include "Function1.H"

//#include "dimensionedScalar.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{
//- Forward declaration of classes
//class electrolyteModels::electrolyteModel;

/*---------------------------------------------------------------------------*\
                           Class varyingDiff Declaration
\*---------------------------------------------------------------------------*/

class varyingDiff
:
    public diffusivityModel
{
    // Private data
      
        //- Index of charge neutralizing ion
        PtrList<label> refSpecieIndex_;

        //- List of concentrations 
        const PtrList<volScalarField>& C_; 

        //- Standard concentration 
        const dimensionedScalar Cstd_; 

        //- List of polynomials for the diffusivity calculation
        PtrList<Function1<scalar>> DVar_; 
  
  
public:

    //- Runtime type information
    TypeName("varyingDiff");


    // Constructors

        //- construct from components
        
        varyingDiff
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte
        );

        varyingDiff
        (
            const dictionary& dict,
            const concReactionThermo& thermo
        );

        //varyingDiff
        //(
        //    const dictionary& dict,
        //    const concReactionThermo& thermo,
        //    const dimensionedScalar& Cstd 
        //);


    // Destructor

        ~varyingDiff()
        {}


    // Member Functions

        //- Calculate diffusivity 
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

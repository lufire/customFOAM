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
    Foam::binaryConductivityModels::linearMixture

Description
    Constant binary conductivity model.

SourceFiles
    linearMixture.C

\*---------------------------------------------------------------------------*/

#ifndef linearMixture_H
#define linearMixture_H

#include "conductivityModel.H"
//#include "dimensionedScalar.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{


/*---------------------------------------------------------------------------*\
                           Class linearMixture Declaration
\*---------------------------------------------------------------------------*/

class linearMixture
:
    public conductivityModel
{
    // Private data
      
        //- Index of charge neutralizing ion
        label balanceIndex_;

        //- List of equivalent conductivites
        const PtrList<label>& z_; 

        //- List of concentrations 
        const PtrList<volScalarField>& C_; 

        //- List of equivalent conductivites
        PtrList<dimensionedScalar> lambda_; 

        //- List of equivalent conductivites
        PtrList<volScalarField> X_; 
  
  
public:

    //- Runtime type information
    TypeName("linearMixture");


    // Constructors

        //- construct from components
        linearMixture
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte 
        );

    // Destructor

        ~linearMixture()
        {}


    // Member Functions
   
        //- Calculate conductivity 
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

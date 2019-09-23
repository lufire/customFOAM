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
    Foam::binaryConductivityModels::speciesConductivity

Description
    Constant binary conductivity model.

SourceFiles
    speciesConductivity.C

\*---------------------------------------------------------------------------*/

#ifndef speciesConductivity_H
#define speciesConductivity_H

#include "transferenceNumberModel.H"
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
                           Class speciesConductivity Declaration
\*---------------------------------------------------------------------------*/

class speciesConductivity
:
    public transferenceNumberModel
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

        //- List of ionic conductivites at infinite solution 
        PtrList<dimensionedScalar> speciesKappa_; 

        //- List of equivalent conductivites
        PtrList<volScalarField> X_; 
  
        //- List of transference numbers 
        //PtrList<volScalarField> t_; 
  
  
public:

    //- Runtime type information
    TypeName("speciesConductivity");


    // Constructors

        //- construct from components
        speciesConductivity
        (
            const dictionary& dict,
            const electrolyteModel& electrolyte 
        );

    // Destructor

        ~speciesConductivity()
        {}


    // Member Functions
   
        // Return transfernce numbers 
        //inline const PtrList<volScalarField>&  t() const
        //{
        //    return t_;
        //}

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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::thixoShearStressSlipVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This velocity inlet boundary condition is applied to patches where the
    pressure is specified.  The inflow velocity is obtained from the flux with
    a direction normal to the patch faces.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            thixoShearStressSlipVelocity;
        phi             phi;
        rho             rho;
        value           uniform 0;
    }
    \endverbatim

Note
    If reverse flow is possible or expected use
    the pressureInletOutletVelocityFvPatchVectorField condition instead.

SeeAlso
    Foam::fixedValueFvPatchField
    Foam::pressureInletOutletVelocityFvPatchVectorField

SourceFiles
    thixoShearStressSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef thixoShearStressSlipVelocityFvPatchVectorField_H
#define thixoShearStressSlipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class thixoShearStressSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class thixoShearStressSlipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    
        //- Type of slip boundary condition 
        word modelName_;
       
        //- Viscosity field name
        word nuName_;
       
        //- Structure field name
        word lambdaName_;
       
        //- First model parameter 
        const scalar alpha_;
        
        //- Second model parameter 
        const scalar beta_;
        
        //- Parameter for thixotropic influence 
        const scalar gamma_;

        //- Density
        dimensionedScalar rho_;

public:

    //- Runtime type information
    TypeName("thixoShearStressSlipVelocity");


    // Constructors
        
        //- Construct from patch and internal field
        thixoShearStressSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thixoShearStressSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given thixoShearStressSlipVelocityFvPatchVectorField
        //  onto a new patch
        thixoShearStressSlipVelocityFvPatchVectorField
        (
            const thixoShearStressSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thixoShearStressSlipVelocityFvPatchVectorField
        (
            const thixoShearStressSlipVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new thixoShearStressSlipVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thixoShearStressSlipVelocityFvPatchVectorField
        (
            const thixoShearStressSlipVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new thixoShearStressSlipVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            ////- Return the viscosity field name
            //const word& nuName() const
            //{
            //    return nuName_;
            //}

            ////- Return reference to the viscosity field name to allow adjustment
            //word& nuName()
            //{
            //    return nuName_;
            //}


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

    //    virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::stickSlipVelocityFvPatchVectorField

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
        type            stickSlipVelocity;
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
    stickSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef stickSlipVelocityFvPatchVectorField_H
#define stickSlipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class stickSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class stickSlipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    
        //- Viscosity field name
        word nuName_;
       
        //- Structure field name
        word lambdaName_;
       
        //-  Slip exponent 
        const scalar alpha_;
        
        //- Linear model parameter 
        const scalar beta_;
        
        //- Model parameter for thixotropic influence
        const scalar gamma_;
        
        //- Critical wall shear stress 
        const scalar sg_;

        //- Constant density
        scalar rho_;

        //- Yield stress of fluid 
        scalar tau0_;

        //- Transition velocity between non-linear and linear slip 
        const scalar Ustar_;

public:

    //- Runtime type information
    TypeName("stickSlipVelocity");


    // Constructors
        
        //- Construct from patch and internal field
        stickSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        stickSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given stickSlipVelocityFvPatchVectorField
        //  onto a new patch
        stickSlipVelocityFvPatchVectorField
        (
            const stickSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        stickSlipVelocityFvPatchVectorField
        (
            const stickSlipVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new stickSlipVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        stickSlipVelocityFvPatchVectorField
        (
            const stickSlipVelocityFvPatchVectorField&,
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
                new stickSlipVelocityFvPatchVectorField(*this, iF)
            );
        }


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

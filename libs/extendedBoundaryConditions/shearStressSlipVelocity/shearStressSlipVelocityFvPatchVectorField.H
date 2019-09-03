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
    Foam::shearStressSlipVelocityFvPatchVectorField

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
        type            shearStressSlipVelocity;
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
    shearStressSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef shearStressSlipVelocityFvPatchVectorField_H
#define shearStressSlipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class shearStressSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class shearStressSlipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    
        //- Viscosity field name
        word nuName_;
       
        //- Slip law factor 
        const scalar factor_;
        
        //- Slip law exponent 
        const scalar exponent_;

        //- Density
        dimensionedScalar rho_;

public:

    //- Runtime type information
    TypeName("shearStressSlipVelocity");


    // Constructors
        
        //- Construct from patch and internal field
        shearStressSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        shearStressSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given shearStressSlipVelocityFvPatchVectorField
        //  onto a new patch
        shearStressSlipVelocityFvPatchVectorField
        (
            const shearStressSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        shearStressSlipVelocityFvPatchVectorField
        (
            const shearStressSlipVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new shearStressSlipVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        shearStressSlipVelocityFvPatchVectorField
        (
            const shearStressSlipVelocityFvPatchVectorField&,
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
                new shearStressSlipVelocityFvPatchVectorField(*this, iF)
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

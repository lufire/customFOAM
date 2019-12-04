/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::linearGradientFvPatchScalarField

Group
    grpInletBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition.

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            linearGradient;
    }
    \endverbatim

SeeAlso
    Foam::fixedGradientFvPatchField

SourceFiles
    linearGradientFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef linearGradientFvPatchScalarFields_H
#define linearGradientFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class linearGradientFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class linearGradientFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Field name of proportional variable 
        const word fieldName_;

public:

    //- Runtime type information
    TypeName("linearGradient");


    // Constructors

        //- Construct from patch and internal field
        linearGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        linearGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given linearGradientFvPatchScalarField onto
        //  a new patch
        linearGradientFvPatchScalarField
        (
            const linearGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        linearGradientFvPatchScalarField
        (
            const linearGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new linearGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        linearGradientFvPatchScalarField
        (
            const linearGradientFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new linearGradientFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the patch pressure gradient field from the given linearCoeffp
        virtual void updateCoeffs(const scalarField& linearCoeffp);

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"

namespace Foam
{
    template<class GradBC>
    inline void setLinearCoeff
    (
        volScalarField::GeometricBoundaryField& bf,
        const FieldField<fvsPatchField, scalar>& linearCoeff
    )
    {
        forAll(bf, patchi)
        {
            if (isA<GradBC>(bf[patchi]))
            {
                refCast<GradBC>(bf[patchi]).updateCoeffs(linearCoeff[patchi]);
            }
        }
    }

    template<class GradBC>
    inline void setLinearCoeff
    (
        volScalarField::GeometricBoundaryField& bf,
        const tmp<FieldField<fvsPatchField, scalar> >& tlinearCoeff
    )
    {
        setLinearCoeff<GradBC>(bf, tlinearCoeff());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

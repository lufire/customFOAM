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
    Foam::currentPotentialConcentrationFvPatchScalarField

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
        type            butlerVolmerCurrent;
        phi             phi;
        rho             rho;
        value           uniform 0;
    }
    \endverbatim

Note
    If reverse flow is possible or expected use
    the pressureInletOutletVelocityFvPatchScalarField condition instead.

SeeAlso
    Foam::fixedValueFvPatchField
    Foam::pressureInletOutletVelocityFvPatchScalarField

SourceFiles
    currentPotentialConcentrationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef currentPotentialConcentrationFvPatchScalarField_H
#define currentPotentialConcentrationFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class currentPotentialConcentrationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class currentPotentialConcentrationFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data
    
        //- Current density field name
        word currentDensityName_;
        
        //- Electric potential field name
        const word phiEName_;
        
        //- Temperature field name
        const word TName_;
        
         //- Species charge number
        label chargeNumber_;
         
        //- Stoichiometric coefficient
        autoPtr<DataEntry<label> > stoichCoeff_;
        
        //- Number of transferred electrons
        autoPtr<DataEntry<label> > nElectrons_;

        //- Rotation matrix to convert vectors to local coordinate system
        //- based on patch normal vector (i.e. local x-coordinate)
        tensorField rotationMatrix_;

        
public:

    //- Runtime type information
    TypeName("currentPotentialConcentration");


    // Constructors
          
        //- Construct from patch and internal field
        currentPotentialConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
        

        //- Construct from patch, internal field and dictionary
        currentPotentialConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given currentPotentialConcentrationFvPatchScalarField
        //  onto a new patch
        currentPotentialConcentrationFvPatchScalarField
        (
            const currentPotentialConcentrationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        currentPotentialConcentrationFvPatchScalarField
        (
            const currentPotentialConcentrationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new currentPotentialConcentrationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        currentPotentialConcentrationFvPatchScalarField
        (
            const currentPotentialConcentrationFvPatchScalarField&,
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
                new currentPotentialConcentrationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& currentDensityName() const
            {
                return currentDensityName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& currentDensityName()
            {
                return currentDensityName_;
            }
            
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        //virtual void operator=(const fvPatchField<scalar>& psf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::fixedCurrentPotentialFvPatchScalarField

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
    fixedCurrentPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedCurrentPotentialFvPatchScalarField_H
#define fixedCurrentPotentialFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class fixedCurrentPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedCurrentPotentialFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data
    
        //- Name of the electrode (Anode or Cathode)
        word electrodeName_;
         
        //- Current density field name
        word jName_;
        
        //- Temperature field name
        word TName_;
        
        //- Oxidant species field name
        word oxidantName_;
        
        //- Reductant species field name
        word reductantName_;
        
        //- Oxidant species refence mass fraction
        scalar CRefOxidant_;
         
         //- Reductant species refence mass fraction
         scalar CRefReductant_;
         
         //- Oxidant species stoichiometric coefficnent
         scalar gammaOxidant_;
         
         //- Reductant species stoichiometric coefficnent
         scalar gammaReductant_;
 
         //- Exchange current density
         scalar jEx_;
         
        //- Equilibrium potential (i.e. open circuit voltage)
        scalar equilibriumPotential_;
        
        //- Electrode potential
        scalar electrodePotential_;        

        //- Anodic charge transfer coefficient
        scalar alphaA_;
        
        //- Cathodic charge transfer coefficient
        scalar alphaC_;

        
public:

    //- Runtime type information
    TypeName("fixedCurrentPotential");


    // Constructors
          
        //- Construct from patch and internal field
        fixedCurrentPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );
        

        //- Construct from patch, internal field and dictionary
        fixedCurrentPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedCurrentPotentialFvPatchScalarField
        //  onto a new patch
        fixedCurrentPotentialFvPatchScalarField
        (
            const fixedCurrentPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedCurrentPotentialFvPatchScalarField
        (
            const fixedCurrentPotentialFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedCurrentPotentialFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedCurrentPotentialFvPatchScalarField
        (
            const fixedCurrentPotentialFvPatchScalarField&,
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
                new fixedCurrentPotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& jName() const
            {
                return jName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& jName()
            {
                return jName_;
            }
            
            //- Return residual of setting Butler-Volmer-Equation equal to zero
            //  according to Newton's method
            //scalar bvNewtonMethodFunction
            //    (const scalar& overP, const label& faceI);
            
            //- Return residual of setting derivative of Butler-Volmer-Equation
            //  equal to zero according to Newton's method
            //scalar bvNewtonMethodDerivative
            //    (const scalar& overP, const label& faceI);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

    //    virtual void operator=(const fvPatchField<scalar>& psf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::butlerVolmerFixedPotentialFvPatchVectorField

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
        type            butlerVolmerFixedPotential;
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
    butlerVolmerFixedPotentialFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef butlerVolmerFixedPotentialFvPatchVectorField_H
#define butlerVolmerFixedPotentialFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class butlerVolmerFixedPotentialFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class butlerVolmerFixedPotentialFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    
        //- Name of the electrode (anode or cathode)
        word electrodeName_;
    
        //- Electric potential field name
        word phiEName_;
        
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
        scalar iEx_;
        
        //- Equilibrium potential (i.e. open circuit voltage)
        scalar eqPotential_;
        
        //- Electrode potential
        scalar electrodePotential_;

        //- Anodic charge transfer coefficient
        scalar alphaA_;
        
        //- Cathodic charge transfer coefficient
        scalar alphaC_;

        //- Number of transferred electrons
        label nElectrons_;

        
public:

    //- Runtime type information
    TypeName("butlerVolmerFixedPotential");


    // Constructors
        
        //- Construct from patch and internal field
        butlerVolmerFixedPotentialFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        butlerVolmerFixedPotentialFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given butlerVolmerFixedPotentialFvPatchVectorField
        //  onto a new patch
        butlerVolmerFixedPotentialFvPatchVectorField
        (
            const butlerVolmerFixedPotentialFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        butlerVolmerFixedPotentialFvPatchVectorField
        (
            const butlerVolmerFixedPotentialFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new butlerVolmerFixedPotentialFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        butlerVolmerFixedPotentialFvPatchVectorField
        (
            const butlerVolmerFixedPotentialFvPatchVectorField&,
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
                new butlerVolmerFixedPotentialFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& phiEName() const
            {
                return phiEName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiEName()
            {
                return phiEName_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef faradayConcentrationFvPatchScalarFields_H
#define faradayConcentrationFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class faradayConcentrationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class faradayConcentrationFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // molar mass ... if we work with mass concentration
    scalar molarMass_;
public:

    //- Runtime type information
    TypeName("faradayConcentration");


    // Constructors

        //- Construct from patch and internal field
        faradayConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        faradayConcentrationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given faradayConcentrationFvPatchScalarField onto
        //  a new patch
        faradayConcentrationFvPatchScalarField
        (
            const faradayConcentrationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        faradayConcentrationFvPatchScalarField
        (
            const faradayConcentrationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new faradayConcentrationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        faradayConcentrationFvPatchScalarField
        (
            const faradayConcentrationFvPatchScalarField&,
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
                new faradayConcentrationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#endif

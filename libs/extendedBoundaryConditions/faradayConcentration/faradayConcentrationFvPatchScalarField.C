#include "faradayConcentrationFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faradayConcentrationFvPatchScalarField::faradayConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    molarMass_(1.0)
{}

// this constructor is called in the solver and decomposePar
Foam::faradayConcentrationFvPatchScalarField::faradayConcentrationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    molarMass_(dict.lookupOrDefault<scalar>("molarMass", 1.0))
{
    if (dict.found("value"))                    // initial value
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}

// probably only needed for mapFields
Foam::faradayConcentrationFvPatchScalarField::faradayConcentrationFvPatchScalarField
(
    const faradayConcentrationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF),
    molarMass_(ptf.molarMass_)
{
    patchType() = ptf.patchType();

    // Map gradient. Set unmapped values and overwrite with mapped ptf
    gradient() = 0.0;
    gradient().map(ptf.gradient(), mapper);

    // Evaluate the value field from the gradient if the internal field is valid
    if (&iF && iF.size())
    {
        scalarField::operator=
        (
            //patchInternalField() + gradient()/patch().deltaCoeffs()
            // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
            // which fails for AMI patches for some mapping operations
            patchInternalField() + gradient()*(patch().nf() & patch().delta())
        );
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value. This
        // constructor is used when reconstructing fields
        this->map(ptf, mapper);
    }
}


Foam::faradayConcentrationFvPatchScalarField::faradayConcentrationFvPatchScalarField
(
    const faradayConcentrationFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    molarMass_(wbppsf.molarMass_)
{}

Foam::faradayConcentrationFvPatchScalarField::faradayConcentrationFvPatchScalarField
(
    const faradayConcentrationFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    molarMass_(wbppsf.molarMass_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faradayConcentrationFvPatchScalarField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    const double F = 96485.3329; // Faraday constant

    const dictionary& fluidDict = this->db().lookupObject<IOdictionary>
    (
                "transportProperties"
    );
    dimensionedScalar z1(fluidDict.lookup("z1"));
    dimensionedScalar D1(fluidDict.lookup("D1"));
    dimensionedScalar z2(fluidDict.lookup("z2"));
// species 1 has to be active ion and 2 has to be passive ion. 
   dimensionedScalar zeff=z1.value()*(z1.value()-z2.value())/z2.value();
 
    // current density
    const fvPatchField<vector>& J = patch().lookupPatchField<volVectorField, scalar>("i");
// species 1: active ion. Species 2: passive ion.
    gradient() = (J / (zeff.value()*F*D1.value()) * molarMass_) & patch().nf();

    fixedGradientFvPatchScalarField::updateCoeffs();
}

void Foam::faradayConcentrationFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("molarMass") << molarMass_ << token::END_STATEMENT << nl;;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        faradayConcentrationFvPatchScalarField
    );
}

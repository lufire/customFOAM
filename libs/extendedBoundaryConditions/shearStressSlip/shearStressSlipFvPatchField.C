/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "shearStressSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::shearStressSlipFvPatchField<Type>::shearStressSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    nuName_("nu"),
    factor_(1.0),
    exponent_(1.0)
{}


template<class Type>
Foam::shearStressSlipFvPatchField<Type>::shearStressSlipFvPatchField
(
    const shearStressSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    nuName_(ptf.nuName_, mapper),
    factor_(ptf.factor_, mapper),
    exponent_(ptf.exponent_, mapper)
{}


template<class Type>
Foam::shearStressSlipFvPatchField<Type>::shearStressSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF),
    nuName_(dict.lookupOrDefault<word>("nuField","nu")),
    factor_(readScalar(dict.lookup("factor"))),
    exponent_(readScalar(dict.lookup("exponent")))
    exponent_(readScalar(dict.lookup("exponent")))
{
    evaluate();
}


template<class Type>
Foam::shearStressSlipFvPatchField<Type>::shearStressSlipFvPatchField
(
    const shearStressSlipFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf),
    nuName_(ptf.nuName_),
    factor_(ptf.factor_),
    exponent_(ptf.exponent_)
{}


template<class Type>
Foam::shearStressSlipFvPatchField<Type>::shearStressSlipFvPatchField
(
    const shearStressSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    nuName_(ptf.nuName_),
    factor_(ptf.factor_),
    exponent_(ptf.exponent_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::shearStressSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchField<Type>::autoMap(m);
    nuName_.autoMap(m);
    factor_.autoMap(m);
    exponent_.autoMap(m);
}


template<class Type>
void Foam::shearStressSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const shearStressSlipFvPatchField<Type>& dmptf =
        refCast<const shearStressSlipFvPatchField<Type> >(ptf);

    nuName_.rmap(dmptf.nuName_, addr);
    factor_.rmap(dmptf.factor_, addr);
    exponent_.rmap(dmptf.exponent_, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::shearStressSlipFvPatchField<Type>::snGrad() const
{
    tmp<Field<Type>> nHat = this->patch().nf();
    const Field<Type> pif(this->patchInternalField());

    const Field<Type>& Uint = this->internalField();
    const GeometricField<Type>& U =
         db().lookupObject<GeometricField<Type> >(Uint.name());

    const GeometricField<Type>& Uprev = U.prevIter();
    const Field<Type>& Uprevp =
        Uprev.boundaryField()[this->patch().name()];

    const tmp<Field<Type> > UprevpI =
          patch().patchInternalField(Uprev.internalField());
    
    const volScalarField& nu = 
        db().lookupObject<volScalarField>(nuName_); 
    const volScalarField& nuPrev = nu.prevIter();
    const scalarField nuPrevp = 
        patch().patchField<volScalarField, scalar>(nuPrev);
    
    const Field<Type> tUprevp = transform(I - sqr(nHat), Uprevp);
    const Field<Type> tUprevpI = transform(I - sqr(nHat), UprevpI);
    //const Field<Type> tUprevp = Uprevp - n*(n & Uprevp);
    //const Field<Type> tUprevpI = UprevpI() - n*(n & UprevpI());
    const Field<Type> nUgradPrev = (tUprevp - tUprevpI)*patch().deltaCoeffs();
    
    const scalarField d = factor_*pow(max(mag(nUgradPrev),SMALL), (exponent_-1.0))
                *pow(nuPrevp, exponent_)*patch().deltaCoeffs();
    //Info << "d: " << d << endl;
    const scalar alpha = 1.0;
    //const Field<Type> Up = (1.0-alpha)*Uprevp + alpha*d/(d+1.0)*tUprevpI;
    const Field<Type> Up = d/(d+1.0)*tUprevpI;
    return
    (
        Up - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::shearStressSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }


    tmp<Field<Type>> nHat = this->patch().nf();

    const Field<Type>& Uint = this->internalField();
    const GeometricField<Type>& U =
         db().lookupObject<GeometricField<Type> >(Uint.name());

    const GeometricField<Type>& Uprev = U.prevIter();
    const Field<Type>& Uprevp =
        Uprev.boundaryField()[this->patch().name()];
    //const fvPatchField<Type>& Uprevp =
    //    this->patch().patchField<Field<Type>, Type>(Uprev);
    //const Field<Type>& Uprev = U.prevIter();

    //const tmp<Field<Type>> UprevpI =
    //    this->patchInternalField();

    const tmp<Field<Type> > UprevpI =
          patch().patchInternalField(Uprev.internalField());
    //const vectorField& UprevpI = tUprevpI();
    //const Field<Type>& UprevpI = tUprevpI();
    
    const volScalarField& nu = 
        db().lookupObject<volScalarField>(nuName_); 
    const volScalarField& nuPrev = nu.prevIter();
    const scalarField nuPrevp = 
        patch().patchField<volScalarField, scalar>(nuPrev);
    
    //tmp<Field<Type>> tUi = this->patchInternalField();
    //tmp<Field<Type>> tn = patch().nf();
    //const Field<Type> n = tn();
    Info << "I: " << I << endl;
    const Field<Type> tUprevp = transform(I - sqr(nHat), Uprevp);
    const Field<Type> tUprevpI = transform(I - sqr(nHat), UprevpI);
    //const Field<Type> tUprevp = Uprevp - n*(n & Uprevp);
    //const Field<Type> tUprevpI = UprevpI() - n*(n & UprevpI());
    const Field<Type> nUgradPrev = (tUprevp - tUprevpI)*patch().deltaCoeffs();
    
    //Info << "factor: " << factor_ << endl;
    //Info << "exponent: " << exponent_ << endl;
    //Info << "deltaCoeffs: " << patch().deltaCoeffs() << endl;
    //Info << "tUprevp: " << tUprevp << endl;
    //Info << "tUprevpI: " << tUprevpI << endl;
    //Info << "nUgradPrev: " << nUgradPrev << endl;
    const scalarField d = factor_*pow(max(mag(nUgradPrev),SMALL), (exponent_-1.0))
                *pow(nuPrevp, exponent_)*patch().deltaCoeffs();
    //Info << "d: " << d << endl;
    const scalar alpha = 1.0;
    //const Field<Type> Up = (1.0-alpha)*Uprevp + alpha*d/(d+1.0)*tUprevpI;
    const Field<Type> Up = d/(d+1.0)*tUprevpI;
    Info << "Patch name: " << patch().name() << endl;
    Info << "tUprevp: " << tUprevp << endl;
    Info << "tUprevpI: " << tUprevpI << endl;
    Info << "Uprevp: " << Uprevp << endl;
    Info << "Up: " << Up << endl;
    Field<Type>::operator=
    (    
        Up
    );

    //fixedValueFvPatchVectorField::updateCoeffs();

    //Field<Type>::operator=
    //(
    //    (1.0 - valueFraction_)
    //   *transform(I - sqr(nHat), this->patchInternalField())
    //);

    transformFvPatchField<Type>::evaluate();
}


//template<class Type>
//Foam::tmp<Foam::Field<Type> >
//Foam::shearStressSlipFvPatchField<Type>::snGradTransformDiag() const
//{
//    const Field<Type> nHat(this->patch().nf());
//    Field<Type> diag(nHat.size());
//
//    diag.replace(vector::X, mag(nHat.component(vector::X)));
//    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
//    diag.replace(vector::Z, mag(nHat.component(vector::Z)));
//
//    return
//        valueFraction_*pTraits<Type>::one
//      + (1.0 - valueFraction_)
//       *transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
//}


template<class Type>
void Foam::shearStressSlipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    nuName_.writeEntry("nu", os);
    factor_.writeEntry("factor", os);
    exponent_.writeEntry("exponent", os);
}


// ************************************************************************* //

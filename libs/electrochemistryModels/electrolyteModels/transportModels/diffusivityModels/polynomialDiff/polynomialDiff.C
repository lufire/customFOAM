/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "polynomialDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "electrolyteModel.H"
//#include "conductivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(polynomialDiff, 0);
addToRunTimeSelectionTable
(
    diffusivityModel,
    polynomialDiff, 
    dictionary
);
addToRunTimeSelectionTable
(
    diffusivityModel,
    polynomialDiff, 
    constant 
);
//addToRunTimeSelectionTable
//(
//    diffusivityModel,
//    polynomialDiff, 
//    Electrolyte 
//);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polynomialDiff::polynomialDiff
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    diffusivityModel(dict, electrolyte),
    refSpecieIndex_(species_.size()),
    C_(electrolyte.thermo().composition().C()),
    Cstd_(electrolyte.Cstd()),
    DPoly_(species_.size())
{
    if(Cstd_.value() <= 0.0) 
    {
        FatalErrorIn
        (
            "diffusivityModel::polynomialDiff::polynomialDiff( "
            "const dictionary&, const electrolyteModel&)"
        )   << "Standard concentration Cstd "
            << Cstd_ << endl 
            << "must be greater than zero." << endl
            << exit(FatalError);
    }
    forAll(species_,specieI)
    {
        refSpecieIndex_.set(specieI, new label(-1));

        dictionary subdict(dict.subDict(species_[specieI]));
        word refSpecie = subdict.lookup("referenceSpecie");

        for(label specieJ = 0; specieJ<species_.size(); ++specieJ)
        {
            if(species_[specieJ] == refSpecie)
            {
                refSpecieIndex_[specieI] = specieJ;
            }
        }
        if(refSpecieIndex_[specieI] == -1)
        {
            FatalErrorIn
            (
                "diffusivityModel::polynomialDiff::polynomialDiff("
                "const dictionary&, const electrolyteModel&)"
            )   << "Unknown reference specie "
                << refSpecie << endl << endl
                << "Valid  reference species are: " << endl
                << species_ 
                << exit(FatalError);
        }
        DPoly_.set
        (
            specieI, 
            new Foam::Function1Types::Polynomial<scalar>
            (
                "diffusivityCoeffs",
                subdict
            )
        );
    }
    update();
}


polynomialDiff::polynomialDiff
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    diffusivityModel(dict, thermo),
    refSpecieIndex_(species_.size()),
    C_(thermo.composition().C()),
    Cstd_(dict.lookup("referenceConcentration")),
    DPoly_(species_.size())
{
    if(Cstd_.value() <= 0.0) 
    {
        FatalErrorIn
        (
            "diffusivityModel::polynomialDiff::polynomialDiff( "
            "const dictionary&, const concReactionThermo&)"
        )   << "Standard concentration Cstd "
            << Cstd_ << endl 
            << "must be greater than zero." << endl
            << exit(FatalError);
    }
    forAll(species_,specieI)
    {
        refSpecieIndex_.set(specieI, new label(-1));

        dictionary subdict(dict.subDict(species_[specieI]));
        word refSpecie = subdict.lookup("referenceSpecie");

        for(label specieJ = 0; specieJ<species_.size(); ++specieJ)
        {
            if(species_[specieJ] == refSpecie)
            {
                refSpecieIndex_[specieI] = specieJ;
            }
        }
        if(refSpecieIndex_[specieI] == -1)
        {
            FatalErrorIn
            (
                "diffusivityModel::polynomialDiff::polynomialDiff("
                "const dictionary&, const concReactionThermo&)"
            )   << "Unknown reference specie "
                << refSpecie << endl << endl
                << "Valid  reference species are: " << endl
                << species_ 
                << exit(FatalError);
        }

        DPoly_.set
        (
            specieI, 
            new Foam::Function1Types::Polynomial<scalar>
            (
                "diffusivityCoeffs",
                subdict
            )
        );
    }
    update();
}


//polynomialDiff::polynomialDiff
//
//(
//    const dictionary& dict,
//    const concReactionThermo& thermo,
//    const dimensionedScalar& Cstd 
//)
//:
//    diffusivityModel(dict, thermo, Cstd),
//    refSpecieIndex_(species_.size()),
//    C_(thermo.composition().C()),
//    Cstd_(Cstd),
//    DPoly_(species_.size())
//{
//    forAll(species_,specieI)
//    {
//        dictionary subdict(dict.subDict(species_[specieI]));
//        word refSpecie = subdict.lookup("referenceSpecie");
//
//        for(label specieJ = 0; specieJ<species_.size(); ++specieJ)
//        {
//            if(species_[specieJ] == refSpecie)
//                refSpecieIndex_[specieI] = specieJ;
//            else
//            {
//                FatalErrorIn
//                (
//                    "diffusivityModel::polynomialDiff::polynomialDiff( "
//                    "const dictionary&, const volScalarField&,"
//                    "const concReactionThermo&, const wordList&)"
//                )   << "Unknown reference specie"
//                    << refSpecie << endl << endl
//                    << "Valid  reference species_ are: " << endl
//                    << species_ 
//                    << exit(FatalError);
//            }
//        }
//        DPoly_.set
//        (
//            specieI, 
//            new Foam::Function1Types::Polynomial<scalar>
//            (
//                "diffusivityCoeffs",
//                subdict
//            )
//        );
//    }
//    update();
//}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void polynomialDiff::update()
{
    forAll(D_,specieI)
    {
        //volScalarField& Di = D_[specieI];
        scalar y = 0.0;
        forAll(D_[specieI], cellI)
        {
            y = C_[refSpecieIndex_[specieI]][cellI]/Cstd_.value();
            D_[specieI][cellI] = DPoly_[specieI].value(y);
        }
        forAll(D_[specieI].boundaryField(), patchI)
        {
            const fvPatch& patch = D_[specieI].boundaryField()[patchI].patch();
            const labelUList& faceCells = patch.faceCells();
            forAll(D_[specieI].boundaryField()[patchI], faceI)
            {
                y = C_[refSpecieIndex_[specieI]].boundaryField()[patchI][faceI]
                    /Cstd_.value();
                scalar Dvalue = DPoly_[specieI].value(y);
                D_[specieI].boundaryField()[patchI][faceI] = Dvalue;
                D_[specieI][faceCells[faceI]] = Dvalue;
            }
        }
    }
}


} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

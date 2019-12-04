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

#include "varyingDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "electrolyteModel.H"
#include "volFields.H"
//#include "conductivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace electrochemistryModels 
{
namespace electrolyteModels
{

defineTypeNameAndDebug(varyingDiff, 0);
addToRunTimeSelectionTable
(
    diffusivityModel,
    varyingDiff, 
    dictionary
);
addToRunTimeSelectionTable
(
    diffusivityModel,
    varyingDiff, 
    constant 
);
//addToRunTimeSelectionTable
//(
//    diffusivityModel,
//    varyingDiff, 
//    Electrolyte 
//);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

varyingDiff::varyingDiff
(
    const dictionary& dict,
    const electrolyteModel& electrolyte 
)
:
    diffusivityModel(dict, electrolyte),
    refSpecieIndex_(species_.size()),
    C_(electrolyte.thermo().composition().C()),
    Cstd_(electrolyte.Cstd()),
    DVar_(species_.size())
{
    if(Cstd_.value() <= 0.0) 
    {
        FatalErrorIn
        (
            "diffusivityModel::varyingDiff::varyingDiff( "
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
                "diffusivityModel::varyingDiff::varyingDiff("
                "const dictionary&, const electrolyteModel&)"
            )   << "Unknown reference specie "
                << refSpecie << endl << endl
                << "Valid  reference species are: " << endl
                << species_ 
                << exit(FatalError);
        }
        DVar_.set
        (
            specieI, 
            Function1<scalar>::New
            (
                "diffusivityCoeffs",
                subdict
            )
        );
    }
    update();
}


varyingDiff::varyingDiff
(
    const dictionary& dict,
    const concReactionThermo& thermo
)
:
    diffusivityModel(dict, thermo),
    refSpecieIndex_(species_.size()),
    C_(thermo.composition().C()),
    Cstd_(dict.lookup("referenceConcentration")),
    DVar_(species_.size())
{
    if(Cstd_.value() <= 0.0) 
    {
        FatalErrorIn
        (
            "diffusivityModel::varyingDiff::varyingDiff( "
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
                "diffusivityModel::varyingDiff::varyingDiff("
                "const dictionary&, const concReactionThermo&)"
            )   << "Unknown reference specie "
                << refSpecie << endl << endl
                << "Valid  reference species are: " << endl
                << species_ 
                << exit(FatalError);
        }

        DVar_.set
        (
            specieI, 
            Function1<scalar>::New
            (
                "diffusivityCoeffs",
                subdict
            )
        );
    }
    update();
}


//varyingDiff::varyingDiff
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
//    DVar_(species_.size())
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
//                    "diffusivityModel::varyingDiff::varyingDiff( "
//                    "const dictionary&, const volScalarField&,"
//                    "const concReactionThermo&, const wordList&)"
//                )   << "Unknown reference specie"
//                    << refSpecie << endl << endl
//                    << "Valid  reference species_ are: " << endl
//                    << species_ 
//                    << exit(FatalError);
//            }
//        }
//        DVar_.set
//        (
//            specieI, 
//            Function1<scalar>::New
//            (
//                "diffusivityCoeffs",
//                subdict
//            )
//        );
//    }
//    update();
//}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void varyingDiff::update()
{
    forAll(D_,specieI)
    {
        //volScalarField& Di = D_[specieI];
        scalar y = 0.0;
        forAll(D_[specieI], cellI)
        {
            y = C_[refSpecieIndex_[specieI]][cellI]/Cstd_.value();
            D_[specieI][cellI] = DVar_[specieI].value(y);
        }
        volScalarField::Boundary& Dboundary = D_[specieI].boundaryFieldRef();
        forAll(Dboundary, patchI)
        {
            const fvPatch& patch = Dboundary[patchI].patch();
            const labelUList& faceCells = patch.faceCells();
            forAll(Dboundary[patchI], faceI)
            {
                y = C_[refSpecieIndex_[specieI]].boundaryField()[patchI][faceI]
                    /Cstd_.value();
                scalar Dvalue = DVar_[specieI].value(y);
                Dboundary[patchI][faceI] = Dvalue;
                D_[specieI][faceCells[faceI]] = Dvalue;
            }
        }
    }
}


} // End namespace electrolyteModels 
} // End namespace electrochemistryModels 
} // End namespace Foam

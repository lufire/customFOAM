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

#include "constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(constant, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            constant, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::constant::constant
(
    const word& name1,
    const word& name2,
    const dictionary& dic,
    const volScalarField& p,
    const volScalarField& T
)
:
    binaryDiffusivityModel(name1, name2, dic, p, T)
{
    dimensionedScalar coefficientValue =
        dimensionedScalar("Dij", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1);

    word couplename = name1 + word("-") + name2;
    if
    (
        dic.subDict("constantCoefficients")
            .lookupEntryPtr(couplename,1,1) == NULL
    )
    {
        couplename = name2 + word("-") + name1;
    }
    
    if(name1 != name2)
    {
        Dvalue_ =
            dimensionedScalar(dic.subDict("constantCoefficients")
                .lookup(couplename)).value();
    }

    /*
    size_t found = name1.find("-");
    if(found == std::string::npos)
    {
        found = name1.find("+");
        if(found == std::string::npos)
        {
            z1_ = 0;
        }
        else
        {
            z1_ = name1.size() - found;
        }
    }
    else
    {
        z1_ = (name1.size() - found)*(-1);
    }
    
    found = name2.find("-");
    if(found == std::string::npos)
    {
        found = name2.find("+");
        if(found == std::string::npos)
        {
            z2_ = 0;
        }
        else
        {
            z2_ = name2.size() - found;
        }
    }
    else
    {
        z2_ = (name2.size() - found)*(-1);
    }
    */
    //Info << "name1: " << name1 << ", z1: " << z1_ << endl;
    //Info << "name2: " << name2 << ", z2: " << z2_ << endl;
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModels::constant::D
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD.ref();

    forAll(T, facei)
    {
        d[facei] = Dvalue_;
    }

    return tD;
}


Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::constant::D() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "D_" + name1_ + "_" + name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -1, 0, 0)
        )
    );

    volScalarField& d = tD.ref();

    forAll(this->T_, celli)
    {
        d[celli] = Dvalue_;
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = Dvalue_;
        }
    }

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

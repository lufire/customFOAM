/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::concMultiComponentMixture

Description
    Foam::concMultiComponentMixture

SourceFiles
    concMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef concMultiComponentMixture_H
#define concMultiComponentMixture_H

#include "concBasicSpecieMixture.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class concMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class concMultiComponentMixture
:
    public concBasicSpecieMixture
{
    // Private Data

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Calculate the mass fractions from the concentrations 
        //  and the mixture density
        void calculateMassFractions();

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Construct as copy (not implemented)
        concMultiComponentMixture(const concMultiComponentMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  mesh and phase name
        concMultiComponentMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&,
            const word&
        );

        //- Construct from dictionary, mesh and phase name
        concMultiComponentMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~concMultiComponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "concMultiComponentMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "concMultiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

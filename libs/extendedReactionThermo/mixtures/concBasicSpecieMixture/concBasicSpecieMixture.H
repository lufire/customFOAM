/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::concBasicSpecieMixture

Description
    Specialization of concBasicMultiComponentMixture for a mixture consisting
    of a number for molecular species.

SourceFiles
    concBasicSpecieMixture.C

\*---------------------------------------------------------------------------*/

#ifndef concBasicSpecieMixture_H
#define concBasicSpecieMixture_H

#include "concBasicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class concBasicSpecieMixture Declaration
\*---------------------------------------------------------------------------*/

class concBasicSpecieMixture
:
    public concBasicMultiComponentMixture
{

public:

    //- Run time type information
    TypeName("concBasicSpecieMixture");

    //- The base class of the mixture
    typedef concBasicSpecieMixture basicMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        concBasicSpecieMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&,
            const word&
        );


    //- Destructor
    virtual ~concBasicSpecieMixture()
    {}


    // Member Functions

        // Per specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar Wi(const label speciei) const = 0;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual scalar Cv
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label speciei) const = 0;

            //- Entropy [J/kg/K]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Helmholtz free energy [J/kg]
            virtual scalar A
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;


        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal diffusivity of enthalpy [kg/m/s]
            virtual scalar alphah
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

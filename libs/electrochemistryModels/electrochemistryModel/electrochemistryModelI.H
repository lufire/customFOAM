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
    ANelectrode WARRANTelectrode; without even the implied warranty of MERCHANTABILITelectrode or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    electrodeou should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

namespace Foam
{
namespace electrochemistryModels 
{

inline const IOdictionary& 
electrochemistryModel::dict() const
{
    return electrochemDict_;
}

inline const electrolyteModels::electrolyteModel& 
electrochemistryModel::electrolyte() const
{
    return electrolyte_();
}

inline electrolyteModels::electrolyteModel& 
electrochemistryModel::electrolyte()
{
    return electrolyte_();
}

inline PtrList<electrodeModels::electrodeModel>&
electrochemistryModel::electrodes()
{
    return electrodes_;
}


inline const PtrList<electrodeModels::electrodeModel>&
electrochemistryModel::electrodes() const
{
    return electrodes_;
}


inline electrodeModels::electrodeModel& 
electrochemistryModel::electrode
(
    const label i
)
{
    return electrodes_[i];
}


inline const electrodeModels::electrodeModel& 
electrochemistryModel::electrode
(
    const label i
) const
{
    return electrodes_[i];
}


inline electrodeModels::electrodeModel& 
electrochemistryModel::electrode
(
    const word& electrodeName
)
{
    return electrodes_[electrodeTable_[electrodeName]];
}


inline const electrodeModels::electrodeModel& 
electrochemistryModel::electrode
(
    const word& electrodeName
) const
{
    return electrodes_[electrodeTable_[electrodeName]];
}

inline bool 
electrochemistryModel::active() 
{
    return active_;
}

} // End namespace electrochemistryModels 
} // End namespace Foam
// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "M.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(M, 0);
    addToRunTimeSelectionTable(functionObject, M, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::M::calc()
{
    if (   foundObject<volScalarField>(k_resolved_name) 
        && foundObject<volScalarField>(k_name) )
    {
        const auto& k_res = lookupObject<volScalarField>(k_resolved_name);
        const auto& k = lookupObject<volScalarField>(k_name);

        return store( resultName_, k_res / (k + k_res) );
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::M::M
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict,
        dict.getOrDefault<word>("k", "k"), "M"),
    k_name(dict.getOrDefault<word>("k", "k")),
    k_resolved_name(dict.getOrDefault<word>("kRes", "resolvedKineticEnergy"))
{
    
}


// ************************************************************************* //

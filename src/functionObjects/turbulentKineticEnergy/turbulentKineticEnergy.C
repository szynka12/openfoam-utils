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

#include "turbulentKineticEnergy.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulentKineticEnergy, 0);
    addToRunTimeSelectionTable(functionObject, turbulentKineticEnergy, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::turbulentKineticEnergy::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        const volVectorField& UMean = lookupObject<volVectorField>("UMean");
        const volVectorField Uprime(U-UMean);
        const volScalarField magUprime(Uprime.component(vector::X)*Uprime.component(vector::X) + 
                                       Uprime.component(vector::Y)*Uprime.component(vector::Y) +
                                       Uprime.component(vector::Z)*Uprime.component(vector::Z));

        return store
        (
            resultName_,
            0.5 * magUprime
        );

    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulentKineticEnergy::turbulentKineticEnergy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName(typeName, "U");
}


// ************************************************************************* //

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

#include "kBudgetCorrelations.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(kBudgetCorrelations, 0);
    addToRunTimeSelectionTable(functionObject, kBudgetCorrelations, dictionary);
    
    word kBudgetCorrelations::reynoldsStressTensorName 
      = "reynoldsStressTensor";
    word kBudgetCorrelations::turbulentDiffusionCorrelationName 
      = "turbulentDiffusionCorrelation";
    word kBudgetCorrelations::velocityPressureCorrelationName 
      = "velocityPressureGradientCorrelation";
    word kBudgetCorrelations::dissipationCorrelationName
      = "dissipationCorrelation";
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::kBudgetCorrelations::calc()
{
    if (   
           foundObject<volVectorField>(fieldName_) 
        && foundObject<volVectorField>(fieldName_ + "Mean")
        && foundObject<volVectorField>(pressureFieldName_ )
        && foundObject<volVectorField>(pressureFieldName_  + "Mean")
       )
    {
        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        const volVectorField& UMean 
          = lookupObject<volVectorField>(fieldName_ + "Mean");
        
        const volScalarField& P = lookupObject<volScalarField>(pressureFieldName_);
        const volScalarField& PMean 
          = lookupObject<volScalarField>(pressureFieldName_ + "Mean");
        
        const tmp<volVectorField> UPrim(U - UMean);
        const tmp<volScalarField> PPrim(P - PMean);
        const tmp<volTensorField> gradUPrim(fvc::grad(UPrim));

        bool result = true;
        result = result 
          && store<volTensorField>(reynoldsStressTensorName, UPrim * UPrim);
        
        result = result 
          && store<volVectorField>(
              turbulentDiffusionCorrelationName, (UPrim & UPrim) * UPrim );

        result = result 
          && store<volScalarField>(
              velocityPressureCorrelationName, UPrim & fvc::grad(PPrim) );

        result = result 
          && store<volScalarField>(
              dissipationCorrelationName, gradUPrim && gradUPrim );

        return result;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::kBudgetCorrelations::kBudgetCorrelations
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, 
        dict.lookupOrDefault<word>("field", "U") ),
    pressureFieldName_(dict.lookupOrDefault<word>("pressure", "p"))
{}


// ************************************************************************* //

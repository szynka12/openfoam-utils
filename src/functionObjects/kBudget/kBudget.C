#include "kBudget.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldAverageItem.H"
#include "Tensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(kBudget, 0);
    addToRunTimeSelectionTable(functionObject, kBudget, dictionary);

    word kBudget::productionName = "production";
    word kBudget::turbulentDissipationName = "turbulentDissipation";
    word kBudget::velocityPressureName  = "velocityPressureGradient";
    word kBudget::molecularDiffusionName  = "molecularDiffusion";
    word kBudget::dissipationName = "dissipation";

}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::kBudget::calc()
{
    if (   
           foundObject<volVectorField>(fieldName_) 
        && foundObject<volVectorField>(fieldName_ + "Mean")
       )
    {
        const volVectorField& UMean 
          = lookupObject<volVectorField>(fieldName_ + "Mean");
        
        const tmp<volTensorField> gradUMean(fvc::grad(UMean));

        const volTensorField reynoldsStressMean
          = lookupObject<volTensorField>(
              kBudgetCorrelations::reynoldsStressTensorName + "Mean");

        bool result = true;
        result = result 
          && store<volScalarField>(productionName, 
              -2 * reynoldsStressMean && gradUMean);

        const volVectorField& turbulentDiffusionCorrelationMean
          = lookupObject<volVectorField>(
                kBudgetCorrelations::turbulentDiffusionCorrelationName + "Mean"
              );

        result = result 
          && store<volScalarField>(turbulentDissipationName, 
              fvc::div(turbulentDiffusionCorrelationMean));

        const volScalarField& velPressGradCor
          = lookupObject<volScalarField>( 
              kBudgetCorrelations::velocityPressureCorrelationName + "Mean" );

        
        result = result 
          && store<volScalarField>(velocityPressureName, -2 * velPressGradCor);

        
        result = result 
          && store<volScalarField>(molecularDiffusionName, 
              nu*fvc::laplacian(tr(reynoldsStressMean) ) );
        
        const volScalarField& dissipationCorrelationMean 
          = lookupObject<volScalarField>(
                kBudgetCorrelations::dissipationCorrelationName + "Mean");

        
        result = result 
          && store<volScalarField>(dissipationName, 
              2 * nu * dissipationCorrelationMean );


        return result;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::kBudget::kBudget
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, 
        dict.lookupOrDefault<word>("field", "U") ),
    nu(0.0)
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    transportProperties.readEntry("nu", nu);
}


// ************************************************************************* //

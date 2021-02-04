/*---------------------------------------------------------------------------*\
Calculates the correlations needed for the TKE budget

\*---------------------------------------------------------------------------*/

#ifndef KBUDGETCORRELATIONS_H
#define KBUDGETCORRELATIONS_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


class kBudgetCorrelations
:
    public fieldExpression
{

    // Private data
        const word pressureFieldName_; 
        

    // Private Member Functions

        //- Calculate the correlation necesarry for calculation fo TKE budget 
        //  field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("kBudgetCorrelations");

    static word reynoldsStressTensorName;      
    static word turbulentDiffusionCorrelationName;      
    static word velocityPressureCorrelationName;      
    static word dissipationCorrelationName;


    // Constructors

        //- Construct from Time and dictionary
        kBudgetCorrelations
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        kBudgetCorrelations(const kBudgetCorrelations&) = delete;

        //- No copy assignment
        void operator=(const kBudgetCorrelations&) = delete;


    //- Destructor
    virtual ~kBudgetCorrelations() = default;


        //- Write the result field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
Calculates the correlations needed for the TKE budget

\*---------------------------------------------------------------------------*/

#ifndef KBUDGET_H
#define KBUDGET_H

#include "fieldExpression.H"
#include "kBudgetCorrelations/kBudgetCorrelations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


class kBudget
:
    public fieldExpression
{

    // Private Member Functions

        //- Calculate the correlation necesarry for calculation fo TKE budget 
        //  field and return true if successful
        virtual bool calc();

        scalar nu;
        


public:

    //- Runtime type information
    TypeName("kBudget");

    static word productionName;
    static word turbulentDissipationName;
    static word velocityPressureName;
    static word molecularDiffusionName;
    static word dissipationName;
    
    // Constructors

        //- Construct from Time and dictionary
        kBudget
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        kBudget(const kBudget&) = delete;

        //- No copy assignment
        void operator=(const kBudget&) = delete;


    //- Destructor
    virtual ~kBudget() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

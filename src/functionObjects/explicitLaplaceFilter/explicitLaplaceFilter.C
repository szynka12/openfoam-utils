#include "explicitLaplaceFilter.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(explicitLaplaceFilter, 0);
    addToRunTimeSelectionTable(functionObject, explicitLaplaceFilter, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::explicitLaplaceFilter::explicitLaplaceFilter
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fields_(dict.get<List<word>>("fields")),
    filter_(mesh_, dict.get<scalar>("widthCoeff"))
{
    read(dict);

    // add fields from the list to the database
    for (const auto& field : fields_)
    {
      addMeanFieldType<scalar>(field);
      addMeanFieldType<vector>(field);
      addMeanFieldType<symmTensor>(field);
      addMeanFieldType<tensor>(field);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::explicitLaplaceFilter::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::explicitLaplaceFilter::execute()
{
    for (const auto& field : fields_)
    {
        filterFieldType<scalar>(field);
        filterFieldType<vector>(field);
        filterFieldType<symmTensor>(field);
        filterFieldType<tensor>(field);
    
    }
    return true;
}


bool Foam::functionObjects::explicitLaplaceFilter::end()
{
    return true;
}


bool Foam::functionObjects::explicitLaplaceFilter::write()
{
    for (const auto& field : fields_)
    {
      const word filtered_field_name = field + "Filtered";
      writeFieldType<scalar>(filtered_field_name);
      writeFieldType<vector>(filtered_field_name);
      writeFieldType<symmTensor>(filtered_field_name);
      writeFieldType<tensor>(filtered_field_name);
    }
    
    return true;
}


// ************************************************************************* //

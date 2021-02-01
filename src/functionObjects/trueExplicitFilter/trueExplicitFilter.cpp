#include "trueExplicitFilter.h"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(trueExplicitFilter, 0);
    addToRunTimeSelectionTable(functionObject, trueExplicitFilter, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::trueExplicitFilter::trueExplicitFilter
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fields_(dict.get<List<word>>("fields")),
    filterList_(mesh_.nCells())
{

    Info<< "Reading field filter\n" << endl;
    volScalarField filterWidth
    (
        IOobject
        (
            "filterWidth",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    if (!filterWidth.headerOk())
    {
        filterWidth 
          = dimensionedScalar("filterWidth", dimless, dict);
    }

    const auto& cell_centers = mesh_.C();
    for (label celli = 0; celli < mesh_.nCells(); celli++)
    {
      filterList_[celli] 
        = filters::FilterBase<filters::SphericalFilter>(
            cell_centers[celli], filterWidth[celli]);

      filterList_[celli].initialise(mesh_);
    }

    read(dict);

    // add fields from the list to the database
    for (const auto& field : fields_)
    {
      addMeanField<scalar>(field);
      addMeanField<vector>(field);
      addMeanField<symmTensor>(field);
      addMeanField<tensor>(field);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::trueExplicitFilter::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::trueExplicitFilter::execute()
{
    for (const auto& field : fields_)
    {
        filterField<scalar>(field);
        filterField<vector>(field);
        filterField<symmTensor>(field);
        filterField<tensor>(field);
    
    }
    return true;
}


bool Foam::functionObjects::trueExplicitFilter::end()
{
    return true;
}


bool Foam::functionObjects::trueExplicitFilter::write()
{
    for (const auto& field : fields_)
    {
      const word filtered_field_name = field + "Filtered";
      
      writeField<scalar>(filtered_field_name);
      writeField<vector>(filtered_field_name);
      writeField<symmTensor>(filtered_field_name);
      writeField<tensor>(filtered_field_name);
    }

    
    return true;
}


// ************************************************************************* //

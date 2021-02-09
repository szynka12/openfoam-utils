#include "cellExplicitFilter.h"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cellExplicitFilter, 0);
    addToRunTimeSelectionTable(functionObject, cellExplicitFilter, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cellExplicitFilter::cellExplicitFilter
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fields_(dict.get<List<word>>("fields")),
    target_time_(
        runTime.rootPath()/runTime.caseName()/runTime.constant(),
        dict.getOrDefault<string>("targetMesh", "targetMesh")),
    mesh_ptr_(nullptr)
{
    
    // read the other mesh here or in the initialisation
    Foam::Info << Foam::nl 
               <<"    cellFilter: Reading target mesh..."  
               << Foam::endl;
    
    mesh_ptr_.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                target_time_.timeName(),
                target_time_,
                Foam::IOobject::MUST_READ
            )
        )
    );
    label n_taget_cells = mesh_ptr_->nCells();
    Foam::Info << "    cellFilter: Target mesh has " 
               << n_taget_cells << " cells." 
               << Foam::nl << Foam::endl;

    filterList_.resize(n_taget_cells);

    using FilterType = filters::UnstructuredMeshFilter<filters::CellFilter>;

    for (label celli = 0; celli < n_taget_cells; celli++)
    {
      filterList_[celli].set_definition(filters::CellFilter(celli, mesh_ptr_)); 

      filterList_[celli].initialise(mesh_);
    }

    // write filter volumes if options selected in the dict
    if (dict.getOrDefault("writeFilterVolume", false))
    {
      Info << "    cellFilter: Computing filter volume..." << endl;
      mesh_ptr_->thisDb().store
      (
          new volScalarField
          (
              IOobject
              (
                  "FilterVolume",
                  mesh_ptr_->thisDb().time().timeName(0),
                  mesh_ptr_->thisDb(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              mesh_ptr_(),
              dimensioned<scalar>(dimVol, Zero)
          )
      );
      
      volScalarField& V
        = mesh_ptr_->thisDb().lookupObjectRef<volScalarField>("FilterVolume");
      
      for(label celli = 0; celli < mesh_ptr_->nCells(); celli++)
      {
          V[celli] = filterList_[celli].V(); 
      }

      V.correctBoundaryConditions();
      
      Info << "    cellFilter: Writing filter volume..." << nl <<endl;;
      
      V.write();
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

bool Foam::functionObjects::cellExplicitFilter::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::cellExplicitFilter::execute()
{
    Foam::Info << Foam::nl 
               <<"    cellFilter: Filtering variables..."  
               << Foam::nl << Foam::endl;

    for (const auto& field : fields_)
    {
        filterField<scalar>(field);
        filterField<vector>(field);
        filterField<symmTensor>(field);
        filterField<tensor>(field);
    }
    return true;
}


bool Foam::functionObjects::cellExplicitFilter::end()
{
    return true;
}


bool Foam::functionObjects::cellExplicitFilter::write()
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

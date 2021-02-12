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
    mesh_ptr_(nullptr),
    divide_by_volume_(true),
    write_volume_field_(dict.getOrDefault("writeFilterVolume", false))
{
    // check the switches
    
    is_parallel_ = Pstream::parRun();

    if (is_parallel_)
    {
      divide_by_volume_ = false;
      write_volume_field_ = true;
    }

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

    
    // Set the definition for the filter
    for (label celli = 0; celli < n_taget_cells; celli++)
    {
      filterList_[celli].set_definition(filters::CellFilter(celli, mesh_ptr_)); 
    }

    
    Info << "    cellFilter: Computing filter volume and stencils..." << endl;
    meshSearch ms(mesh_);
    for (label celli = 0; celli < n_taget_cells; celli++)
    {
      filterList_[celli].initialise(ms);
    }
    
    Info << "    cellFilter: Done!" << endl;
    
    // If we are parallel then volume will be saved in constant so we will know
    // that it is neccesary for reconstruction
    fileName volume_path;
    if (is_parallel_)
    {
      volume_path = mesh_ptr_->thisDb().time().constant();
    }
    else
    {
      volume_path = mesh_ptr_->thisDb().time().timeName(0);
    }

    
    auto volume_field_header = IOobject("FilterVolume", volume_path,
                mesh_ptr_->thisDb(), IOobject::NO_READ, IOobject::NO_WRITE);

    mesh_ptr_->thisDb().store(
        new volScalarField(
          volume_field_header,mesh_ptr_(), dimensioned<scalar>(dimless, Zero))
    );
    
    // Fill the volume field
    volScalarField& V
      = mesh_ptr_->thisDb().lookupObjectRef<volScalarField>("FilterVolume");
    
    for(label celli = 0; celli < mesh_ptr_->nCells(); celli++)
    {
        V[celli] = filterList_[celli].V(); 
    }

    V.correctBoundaryConditions();
    
    
    // write filter volumes if options selected in the dict
    if (write_volume_field_ || is_parallel_)
    {
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

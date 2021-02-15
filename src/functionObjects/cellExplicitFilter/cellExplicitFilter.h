#ifndef CELLEXPLICITFILTER_H
#define CELLEXPLICITFILTER_H

#include "fvMeshFunctionObject.H"
#include "fvPatchField.H"
#include "GeometricField.H"
#include "Time.H"
#include "UnstructuredMeshFilter.h"
#include "CellFilter.h"
#include "Pstream.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class explicitLaplaceFilter Declaration
\*---------------------------------------------------------------------------*/

class cellExplicitFilter 
:
    public fvMeshFunctionObject
{
    // Private Data
    
        using FilterType = filters::UnstructuredMeshFilter<filters::CellFilter>;
    
        // List of fields that needs to be filtered
        const List<word> fields_;
        
        // Time of the target mesh
        mutable Foam::Time target_time_;
        
        // Target mesh
        Foam::autoPtr<Foam::fvMesh> mesh_ptr_;
        
        // List of filters for each cell of target mesh
        List<FilterType> filterList_;

        // Switches
        bool divide_by_volume_;
        bool write_volume_field_;
        bool is_parallel_;
        bool initialised_;
        bool cache_serial_;
        
    // Helper functions
        void cache_serial() const;

        //- Add mean average field to database
        template<class Type>
        void addMeanField(const word& fieldname);

        //- Write field
        template<class Type>
        void writeField( const word& fieldName ) const;
        
        //- Filter field 
        template<class Type>
        void filterField( const word& fieldName ) const;
public:

    //- Runtime type information
    TypeName("cellExplicitFilter");


    // Constructors

        //- Construct from Time and dictionary
        cellExplicitFilter
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        cellExplicitFilter(const cellExplicitFilter&) = delete;

        //- No copy assignment
        void operator=(const cellExplicitFilter&) = delete;


    //- Destructor
    virtual ~cellExplicitFilter() = default;


    // Member Functions

        //- Read the explicitLaplaceFilter data
        virtual bool read(const dictionary& dict);
        
        virtual void calculate_weights();

        virtual void initialise();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the explicitLaplaceFilter
        virtual bool write();
};

template<class Type>
void cellExplicitFilter::addMeanField(const word& field_name)
{
  typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
  
  if (!foundObject<VolFieldType>(field_name)) { return; }
  
  word filtered_field_name = field_name + "Filtered";
  
  if (mesh_ptr_->foundObject<VolFieldType>(filtered_field_name)) {}
  else if (mesh_ptr_->thisDb().found(filtered_field_name))
  {
      FatalError << "    Cannot allocate average field " 
          << filtered_field_name
          << " since an object with that name already exists."
          << endl;
  }
  else
  {
      const VolFieldType& base_field = lookupObject<VolFieldType>(field_name);

      // Store on registry
      mesh_ptr_->thisDb().store
      (
          new VolFieldType
          (
              IOobject
              (
                  filtered_field_name,
                  mesh_ptr_->thisDb().time().timeName(
                    obr().time().startTime().value()),
                  mesh_ptr_->thisDb(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              mesh_ptr_(),
              dimensioned<Type>(base_field.dimensions(), Zero)
          )
      );
    
  }

  if (mesh_ptr_->foundObject<VolFieldType>(filtered_field_name))
  {
    Foam::Info << "    Field " 
               << filtered_field_name << " added succesfully!" << endl;
  }
  
}

template<class Type>
void cellExplicitFilter::writeField( const word& fieldName ) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    target_time_.setTime(time());
      
    if (mesh_ptr_->foundObject<VolFieldType>(fieldName))
    {
        VolFieldType& f = mesh_ptr_->lookupObjectRef<VolFieldType>(fieldName);
        
        //if (divide_by_volume_)
        //{
          //const volScalarField& V 
            //= mesh_ptr_->lookupObject<volScalarField>("FilterVolume");
          
          //for(label celli = 0; celli < mesh_ptr_->nCells(); celli++)
          //{            
            //f[celli] =  f[celli] * ( 1.0 / (V[celli] + SMALL) );
          //}
        //}

        f.write();
    }
    
    
}

template<class Type>
void cellExplicitFilter::filterField( const word& fieldName ) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    
    if (foundObject<VolFieldType>(fieldName))
    {
      
      Info << "                Filtering field " << fieldName << endl;
      const VolFieldType& f_original = lookupObject<VolFieldType>(fieldName);
      
      VolFieldType& f_filtered 
        = mesh_ptr_->thisDb().lookupObjectRef<VolFieldType>(fieldName + "Filtered");
      
      for(label celli = 0; celli < mesh_ptr_->nCells(); celli++)
      {
          f_filtered[celli] 
            = filterList_[celli].localConvolution<Type>(
                f_original, divide_by_volume_);
      }
      
      f_filtered.correctBoundaryConditions();
      
    }
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

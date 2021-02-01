#ifndef SPHERICALEXPLICITFILTER_H
#define SPHERICALEXPLICITFILTER_H

#include "fvMeshFunctionObject.H"
#include "fvPatchField.H"
#include "GeometricField.H"
#include "SphericalFilter.h"
#include "FilterBase.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class explicitLaplaceFilter Declaration
\*---------------------------------------------------------------------------*/

class trueExplicitFilter 
:
    public fvMeshFunctionObject
{
    // Private Data
        const List<word> fields_;
        
      
        List<filters::FilterBase<filters::SphericalFilter>> filterList_;
    // Helper functions

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
    TypeName("trueExplicitFilter");


    // Constructors

        //- Construct from Time and dictionary
        trueExplicitFilter
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        trueExplicitFilter(const trueExplicitFilter&) = delete;

        //- No copy assignment
        void operator=(const trueExplicitFilter&) = delete;


    //- Destructor
    virtual ~trueExplicitFilter() = default;


    // Member Functions

        //- Read the explicitLaplaceFilter data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the explicitLaplaceFilter
        virtual bool write();
};

template<class Type>
void trueExplicitFilter::addMeanField(const word& field_name)
{
  typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

  if (!foundObject<VolFieldType>(field_name)) { return; }
  
  word filtered_field_name = field_name + "Filtered";
  
  if (foundObject<VolFieldType>(filtered_field_name)) {}
  else if (obr().found(filtered_field_name))
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
      obr().store
      (
          new VolFieldType
          (
              IOobject
              (
                  filtered_field_name,
                  obr().time().timeName(obr().time().startTime().value()),
                  obr(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
              ),
              1.0*base_field
          )
      );

  }
}

template<class Type>
void trueExplicitFilter::writeField( const word& fieldName ) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (foundObject<VolFieldType>(fieldName))
    {
        const VolFieldType& f = lookupObject<VolFieldType>(fieldName);
        f.write();
    }
    
}

template<class Type>
void trueExplicitFilter::filterField( const word& fieldName ) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    
    if (foundObject<VolFieldType>(fieldName))
    {

      const VolFieldType& f_original = lookupObject<VolFieldType>(fieldName);
      
      VolFieldType& f_filtered 
        = obr().lookupObjectRef<VolFieldType>(fieldName + "Filtered");
      
      for(label celli = 0; celli < mesh_.nCells(); celli++)
      {
          f_filtered[celli] 
            = filterList_[celli].localConvolution<Type>(f_original);
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

#ifndef explicitLaplaceFilter_H
#define explicitLaplaceFilter_H

#include "fvMeshFunctionObject.H"
#include "LESfilter.H"
#include "laplaceFilter.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class explicitLaplaceFilter Declaration
\*---------------------------------------------------------------------------*/

class explicitLaplaceFilter
:
    public fvMeshFunctionObject
{
    // Private Data

    List<word> fields_;
    laplaceFilter filter_;

    // Helper functions

        //- Add mean average field to database
        template<class Type>
        void addMeanFieldType(const word& fieldname);

        //- Write field
        template<class Type>
        void writeFieldType( const word& fieldName ) const;
        
        //- Filter field 
        template<class Type>
        void filterFieldType( const word& fieldName ) const;
public:

    //- Runtime type information
    TypeName("explicitLaplaceFilter");


    // Constructors

        //- Construct from Time and dictionary
        explicitLaplaceFilter
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        explicitLaplaceFilter(const explicitLaplaceFilter&) = delete;

        //- No copy assignment
        void operator=(const explicitLaplaceFilter&) = delete;


    //- Destructor
    virtual ~explicitLaplaceFilter() = default;


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
void explicitLaplaceFilter::addMeanFieldType(const word& field_name)
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
              1*base_field
          )
      );

  }
}

template<class Type>
void explicitLaplaceFilter::writeFieldType( const word& fieldName ) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (foundObject<VolFieldType>(fieldName))
    {
        const VolFieldType& f = lookupObject<VolFieldType>(fieldName);
        f.write();
    }
    
}

template<class Type>
void explicitLaplaceFilter::filterFieldType( const word& fieldName ) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    
    if (foundObject<VolFieldType>(fieldName))
    {

      const VolFieldType& f_original = lookupObject<VolFieldType>(fieldName);
      tmp<VolFieldType> t_f_original(f_original);
      VolFieldType& f_filtered 
        = obr().lookupObjectRef<VolFieldType>(fieldName + "Filtered");

      f_filtered = filter_(t_f_original);

      f_filtered.correctBoundaryConditions();
    
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

#ifndef FILTERBASE_H
#define FILTERBASE_H

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "GeometricField.H"
#include "fvMesh.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace filters
{
  template <class FilterDefinition>
    class FilterBase
    : 
      public FilterDefinition
  {
    protected:

    // Private properties
    const fvMesh& mesh_;

    // Cells that will take part in filtering
    List<label> cells_;

    // Weights corresponding to each cell. w_i = V_i * v_i * G_i
    // where G is the value of the filter,  V is the cell volume and v is the
    // optional (usually equal to 1) modifier
    List<scalar> weights_;
    
    // Compute the weight for the cell 
    //virtual scalar weight(const label cell_label) = 0;
    virtual scalar weight(const label cell_label)
    {
      scalar g = G(mesh_.cellCentres()[cell_label]);
      
      if (g == 0.0) {return 0;}
      else
      {
        return g * mesh_.V()[cell_label];
      }
    }


    template<class Type>
    Type filterFieldType( 
        const GeometricField<Type, fvPatchField, volMesh>& field ) const
    {
      Type value(Zero);
      forAll(weights_, i) { value += weights_[i] * field[i]; }
      return value;
    }

    public:

    FilterBase(const fvMesh& mesh, const vector& center, const scalar& width);
    
    void initialise()
    {
      for (label celli = 0; celli < mesh_.nCells(); celli++)
      {
        scalar w = weight(celli);
        if (w > 0.0)
        {
          weights_.append(w);
          cells_.append(celli);
        }
      }
    }

  };

  template<class FilterDefinition>
  FilterBase<FilterDefinition>::FilterBase(
      const fvMesh& mesh, const vector& center, const scalar& width)
  :
    FilterDefinition(center, width),
    mesh_(mesh)
  {
    initialise();
  }


} //namespace filters
} //namespace Foam
} //namespace functionObjectsn

#endif

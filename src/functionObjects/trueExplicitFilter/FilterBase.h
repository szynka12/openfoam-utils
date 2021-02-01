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

    // Cells that will take part in filtering
    List<label> cells_;

    // Weights corresponding to each cell. w_i = V_i * v_i * G_i
    // where G is the value of the filter,  V is the cell volume and v is the
    // optional (usually equal to 1) modifier
    List<scalar> weights_;

    scalar volume_;
    
    //Compute the weight for the cell
    //virtual scalar weight(const label cell_label) = 0;
    virtual std::pair<bool, scalar> 
      weight(const label cell_label, const fvMesh& mesh)
    {
      scalar g = FilterBase::G(mesh.cellCentres()[cell_label]);
      
      if (g < SMALL) {return std::make_pair(false, g);}
      else
      {
        return std::make_pair(true, g * mesh.cellVolumes()[cell_label]);
      }
    }



    public:
    FilterBase() : FilterDefinition(Zero, Zero) {};
    FilterBase(const vector& center, const scalar& width);
    
    void initialise(const fvMesh& mesh)
    {
      for (label celli = 0; celli < mesh.nCells(); celli++)
      {
        std::pair<bool, scalar> w = weight(celli, mesh);
        if (w.first)
        {
          weights_.append(w.second);
          cells_.append(celli);

          volume_+=mesh.cellVolumes()[celli];
        }
      }
    }
    
    template<class Type>
    Type localConvolution( 
        const GeometricField<Type, fvPatchField, volMesh>& field ) const
    {
      Type value(Zero);
      forAll(weights_, i) { value += weights_[i] * field[cells_[i]]; }
      return value/volume_;
    }

  };

  template<class FilterDefinition>
  FilterBase<FilterDefinition>::FilterBase(
      const vector& center, const scalar& width)
  :
    FilterDefinition(center, width),
    volume_(0.0)
  { }


} //namespace filters
} //namespace Foam
} //namespace functionObjectsn

#endif

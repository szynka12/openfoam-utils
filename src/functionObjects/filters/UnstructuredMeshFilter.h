#ifndef UNSTRUCTUREDMESHFILTER_H
#define UNSTRUCTUREDMESHFILTER_H

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
    class UnstructuredMeshFilter
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


    public:

    UnstructuredMeshFilter() : FilterDefinition(0, autoPtr<fvMesh>(nullptr)) {};
    
    UnstructuredMeshFilter(label cell, const autoPtr<fvMesh>& mesh) 
      : FilterDefinition(cell, mesh) 
    {};
    
    // copy assignment
    virtual UnstructuredMeshFilter<FilterDefinition>& operator=(
        const UnstructuredMeshFilter<FilterDefinition>& other)
    {
        Info << "This is happening" << endl;  

        // Guard self assignment
        if (this == &other)
            return *this;
        
        cells_ = other.cells_;
        weights_ = other.weights_;
        volume_ = other.volume_;

        return *this;
    }

    scalar V() const
    {
      return volume_;
    }

    void set_definition(FilterDefinition def)
    {
      FilterDefinition::operator=(def);
    }

    void initialise(const fvMesh& mesh)
    {
      for (label celli = 0; celli < mesh.nCells(); celli++)
      {
        const vector& evaluation_point = mesh.cellCentres()[celli];
        const scalar V = mesh.cellVolumes()[celli];

        if (FilterDefinition::in_range(evaluation_point))
        {
          weights_.append(
              FilterDefinition::G(evaluation_point) * V);
          
          cells_.append(celli);

          volume_ += V;
        }
      }
    }

    template<class Type>
    Type localConvolution( 
        const GeometricField<Type, fvPatchField, volMesh>& field, 
        const bool divide_by_volume = true ) const
    {
      Type value(Zero);
      forAll(weights_, i) { value += weights_[i] * field[cells_[i]]; }
      
      if (divide_by_volume)
        value /= (volume_ + SMALL);
      
      return value;
    }
    
    template<class Type>
    Type localConvolution(
        const Type uniform_value, const bool divide_by_volume = true) const
    {
      Type value(Zero);
      forAll(weights_, i) { value += weights_[i] * uniform_value; }
      
      if (divide_by_volume)
        value /= (volume_ + SMALL);
      
      return value;
    }
    

  };



} //namespace filters
} //namespace Foam
} //namespace functionObjectsn

#endif

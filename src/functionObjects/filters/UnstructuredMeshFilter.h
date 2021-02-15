#ifndef UNSTRUCTUREDMESHFILTER_H
#define UNSTRUCTUREDMESHFILTER_H

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "GeometricField.H"
#include "fvMesh.H"
#include "List.H"
#include "FIFOStack.H"
#include "meshSearch.H"

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

    void add_cell(const label cell, const vector point, const scalar volume);

    void add_to_queue(
        const label cell, FIFOStack<label>& stack, boolList& marked);

    void add_cell_nei_to_queue(const label cell, const labelListList& cell_nei,
        FIFOStack<label>& stack, boolList& marked);

    public:

    UnstructuredMeshFilter() 
      : 
    FilterDefinition(0, autoPtr<fvMesh>(nullptr)) {};
    
    UnstructuredMeshFilter(label cell, const autoPtr<fvMesh>& mesh) 
      : 
    FilterDefinition(cell, mesh) 
    {};
    
    // copy assignment
    //virtual UnstructuredMeshFilter<FilterDefinition>& operator=(
        //const UnstructuredMeshFilter<FilterDefinition>& other)
    //{
        //Info << "This is happening" << endl;  

        //// Guard self assignment
        //if (this == &other)
            //return *this;
        
        //cells_ = other.cells_;
        //weights_ = other.weights_;
        //volume_ = other.volume_;

        //return *this;
    //}

    inline scalar V() const { return volume_; }

    inline List<scalar> weights() const {return weights_;};
    inline List<label> cells() const {return cells_;};

    void set_definition(FilterDefinition def);

    void initialise(const meshSearch& ms);
    void initialise(const List<scalar>& weights, const List<label>& cells, 
        const fvMesh& mesh);


    template<class Type>
    Type localConvolution( 
        const GeometricField<Type, fvPatchField, volMesh>& field, 
        const bool divide_by_volume = true ) const
    {
      Type value(Zero);
      forAll(weights_, i) { value += weights_[i] * field[cells_[i]]; }
      
      if (divide_by_volume)
        value /= (volume_ + VSMALL);
      
      return value;
    }
    
    //template<class Type>
    //Type localConvolution(
        //const Type uniform_value, const bool divide_by_volume = true) const
    //{
      //Type value(Zero);
      //forAll(weights_, i) { value += weights_[i] * uniform_value; }
      
      //if (divide_by_volume)
        //value /= (volume_ + VSMALL);
      
      //return value;
    //}
    

  };

template <class FilterDefinition>
void UnstructuredMeshFilter<FilterDefinition>::add_cell(
    const label cell, const vector point, const scalar volume)
{
  weights_.append(FilterDefinition::G(point) * volume);
  cells_.append(cell);
  volume_ += volume;
};

template <class FilterDefinition>
void UnstructuredMeshFilter<FilterDefinition>::add_to_queue(
    const label cell, FIFOStack<label>& stack, boolList& marked)
{
  if (!marked[cell])
  {
    marked[cell] = true;
    stack.push(cell);
  }
}


template <class FilterDefinition>
void UnstructuredMeshFilter<FilterDefinition>::add_cell_nei_to_queue(
    const label cell, const labelListList& cell_nei, 
    FIFOStack<label>& stack, boolList& marked)
{
  forAll(cell_nei[cell], ni)
  {
    add_to_queue(cell_nei[cell][ni], stack, marked);
  }
}


template <class FilterDefinition>
void UnstructuredMeshFilter<FilterDefinition>::set_definition(FilterDefinition def)
{
  FilterDefinition::operator=(def);
}

template <class FilterDefinition>
void UnstructuredMeshFilter<FilterDefinition>::initialise(
    const List<scalar>& weights, const List<label>& cells, const fvMesh& mesh)
{
  weights_ = weights;
  cells_ = cells;

  const scalarField& cell_volumes = mesh.cellVolumes();
  
  forAll(cells_, i)
  {
    volume_ += cell_volumes[cells_[i]];
  }
}

template <class FilterDefinition>
void UnstructuredMeshFilter<FilterDefinition>::initialise(const meshSearch& ms)
{
  const polyMesh& mesh = ms.mesh();

  // Check if the mesh is even overlaping with the filter. We need two
  // bounding boxes for that and ceck for overlap. If it is not there, then
  // we have no need to check further. This should only make parallel run
  // faster and serial run slightly slower.
  if (!FilterDefinition::bounding_box().overlaps(mesh.bounds()))
    return;
  
  const vectorField& cell_centres = mesh.cellCentres();
  const scalarField& cell_volumes = mesh.cellVolumes();

  const auto& cell_nei = mesh.cellCells();

  // Loop over the cells to find a first cell (or do a tree search) that is 
  // actually in range.
  // Than insert it's neighbours into the FIFO stack. With that we loop over
  // the stack inserting the neighbours, until there are no left.
  FIFOStack<label> cells_queue;
  boolList queued_cells(mesh.nCells(), false);
  
  /*
  for (label celli = 0; celli < mesh.nCells(); celli++)
  {
    const vector& evaluation_point = cell_centres[celli];

    if (FilterDefinition::in_range(evaluation_point))
    {
      add_cell(celli, evaluation_point, cell_volumes[celli]);
      
      //mark this cell as queued
      queued_cells[celli] = true;
      add_cell_nei_to_queue(celli, cell_nei, cells_queue, queued_cells);
      break;
    }
  }
  */

  // Alternatively we use findCell function
  /*
  label celli = mesh.findCell(FilterDefinition::center());
  
  if (celli == -1 )
    return;
  
  //mark this cell as queued
  queued_cells[celli] = true;
  add_cell(celli, cell_centres[celli], cell_volumes[celli]);
  add_cell_nei_to_queue(celli, cell_nei, cells_queue, queued_cells);
  */
  
  // So, I have no idea how that will scale, it seems that it is bullshit
  // slow in parallel cases, and looping over the whole mesh is just faster.
  // It might be due to communication but I am unsure why. Let's see how it
  // will scale.
  //meshSearch ms(mesh);
  
  label celli = ms.findNearestCell(FilterDefinition::center());
  
  if (celli == -1 )
    return;
  
  //mark this cell as queued
  queued_cells[celli] = true;
  add_cell(celli, cell_centres[celli], cell_volumes[celli]);
  add_cell_nei_to_queue(celli, cell_nei, cells_queue, queued_cells);
  
  while(cells_queue.size() != 0)
  {
    label celli = cells_queue.pop();

    const vector& evaluation_point = cell_centres[celli];

    if (FilterDefinition::in_range(evaluation_point))
    {
      add_cell(celli, evaluation_point, cell_volumes[celli]);
      add_cell_nei_to_queue(celli, cell_nei, cells_queue, queued_cells);
    }
  }
}

} //namespace filters
} //namespace Foam
} //namespace functionObjectsn

#endif

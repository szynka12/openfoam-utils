#ifndef CELLFILTER_H
#define CELLFILTER_H

#include "vector.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"
#include "tmp.H"

namespace Foam
{
namespace functionObjects
{
namespace filters
{

class CellFilter
{
  //Private
    label cell_;
    const fvMesh* target_mesh_;

  // Public
    public:
      CellFilter(const label& cell, const autoPtr<fvMesh>& mesh)
      :
        cell_(cell),
        target_mesh_(mesh.get())
      {
        if (target_mesh_ == nullptr)
        {
          FatalError << "Filter for cell " << cell_ << " has nullptr!" << endl;
        }
      }
     
      virtual ~CellFilter(){};

      virtual inline bool in_range(const vector& evaluation_point) const
      {
        if (!target_mesh_) 
          FatalError << "Filter for cell " << cell_ << " has nullptr!" << endl;
       
        return target_mesh_->pointInCell(evaluation_point, cell_);
      }

      virtual inline boundBox bounding_box() const 
      {
        auto cell = target_mesh_->cellPoints(cell_);
        List< point > points(cell.size());
        forAll(cell, celli) { points[celli] = target_mesh_->points()[celli]; }
        return boundBox(points);
      }

      virtual inline scalar G(const vector& evaluation_point) const
      {
        return 1;
      }

      virtual vector center()
      {
        return target_mesh_->cellCentres()[cell_];
      }

      // copy assignment
      virtual CellFilter& operator=(const CellFilter& other)
      {
          // Guard self assignment
          if (this == &other)
              return *this;
          
          cell_ = other.cell_;

          // I don 't own this pointer so nothing has to be done with it;
          target_mesh_ = other.target_mesh_;

          return *this;
      }
  };

}
}
}


#endif

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
      CellFilter(const label& cell, const autoPtr<fvMesh>& mesh);
     
      virtual ~CellFilter(){};

      virtual bool in_range(const vector& evaluation_point) const;

      virtual boundBox bounding_box() const;

      virtual inline scalar G(const vector& evaluation_point) const
      {
        return 1;
      }

      virtual vector center()
      {
        return target_mesh_->cellCentres()[cell_];
      }

      // copy assignment
      virtual CellFilter& operator=(const CellFilter& other);
  };

}
}
}


#endif

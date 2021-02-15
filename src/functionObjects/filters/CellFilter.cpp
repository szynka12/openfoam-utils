#include "CellFilter.h"

namespace Foam
{
namespace functionObjects
{
namespace filters
{

  CellFilter::CellFilter(const label& cell, const autoPtr<fvMesh>& mesh)
   :
     cell_(cell),
     target_mesh_(mesh.get())
   {
     if (target_mesh_ == nullptr)
     {
       FatalError << "from CellFilter: Filter for cell " 
                  << cell_ << " has nullptr!" << endl;
     }
   }


   bool CellFilter::in_range(const vector& evaluation_point) const
   {
     if (!target_mesh_) 
       FatalError << "from in_range: Filter for cell " << cell_ << " has nullptr!" << endl;
    
     return target_mesh_->pointInCell(evaluation_point, cell_);
   }


   boundBox CellFilter::bounding_box() const 
   {
     auto cell = target_mesh_->cellPoints(cell_);
     List< point > points(cell.size());
     forAll(cell, celli) { points[celli] = target_mesh_->points()[celli]; }
     return boundBox(points);
   }
  

   CellFilter& CellFilter::operator=(const CellFilter& other)
   {
       // Guard self assignment
       if (this == &other)
           return *this;
       
       cell_ = other.cell_;

       // I don 't own this pointer so nothing has to be done with it;
       target_mesh_ = other.target_mesh_;

       return *this;
   }

}
}
}


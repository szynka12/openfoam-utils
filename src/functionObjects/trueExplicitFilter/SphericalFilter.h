#ifndef SPHERICALFILTER_H
#define SPHERICALFILTER_H

#include "FilterDefinitionBase.h"
#include "mathematicalConstants.H"

namespace Foam
{
namespace functionObjects
{
namespace filters
{

class SphericalFilter
  :
  public FilterDefinitionBase
{
  
  // Public
    public:
      SphericalFilter(const vector& center, const scalar& width)
        :
        FilterDefinitionBase(center, width)

      {}
     
      virtual ~SphericalFilter(){};

      virtual scalar G(const vector& evaluation_point) const
      {
        return mag(evaluation_point - center_) < (width_ / 2)
          ? 1
          : 0;
      }

};

}
}
}


#endif

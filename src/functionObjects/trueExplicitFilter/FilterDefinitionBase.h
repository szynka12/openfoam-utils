#ifndef FILTERDEFINITIONBASE_H
#define FILTERDEFINITIONBASE_H

#include "scalar.H"
#include "vector.H"

namespace Foam
{
namespace functionObjects
{
namespace filters
{

class FilterDefinitionBase
{
  
  //Protected
    protected:
      // Center of the filter
      vector center_;
      
      // The width of the filter 
      scalar width_;

  // Public
    public:
      FilterDefinitionBase(const vector& center, const scalar& width)
        :
      center_(center),
      width_(width)
      {};

      virtual ~FilterDefinitionBase(){};
      
      virtual scalar G(const vector& evaluation_point) const = 0;

};

}
}
}


#endif

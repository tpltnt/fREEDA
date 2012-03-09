// This is the parent class for all analysis classes
// by Carlos E. Christoffersen

#ifndef Analysis_h
#define Analysis_h 1

#include "../network/Circuit.h"
#include "../containers.h"
#include "Amesos.h"

class Analysis : public NetListItem
{
	public:
  Analysis(ItemInfo* ainfo, ParmInfo* param_desc, const int& numparms);

  virtual ~Analysis();

  virtual void run(Circuit* cir) = 0;
};

#endif


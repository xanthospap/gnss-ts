#ifndef __GEODESY__
#define __GEODESY__

#include "geoconst.hpp"

namespace ngpt {

/// Given a vector in (n, e, u), compute the distance, azimouth and zenith 
/// distance.
void top2daz(double,  double,  double, double&, double&, double&);

} // end namespace geodesy

#endif

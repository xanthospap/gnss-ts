#ifndef _GEODESY_
#define _GEODESY_

#include "geoconst.hpp"

namespace geodesy {

void top2daz(const double& north,const double& east,const double& up,
    double& distance, double& azimouth, double& zenith);

} // end namespace geodesy

#endif

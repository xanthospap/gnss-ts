#ifndef __GEODESY__
#define __GEODESY__

#include "geoconst.hpp"

namespace ngpt
{

/// Given a vector in (n, e, u), compute the distance, azimouth and zenith 
/// distance.
void top2daz(double,  double,  double, double&, double&, double&);

namespace detail
{
void
car2top_matrix(double sinf, double sinl, double cosf, double cosl, double* coef)
noexcept;

void
car2top_cov_matrix(double sin2f, double sin2l, double cos2f, double cos2l, double* coef)
noexcept;
} // end namespace detail

} // end namespace geodesy

#endif

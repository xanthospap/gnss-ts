#ifndef __NGPT_GEO_DETAILS__
#define __NGPT_GEO_DETAILS__

#include <cmath>

namespace ngpt
{
namespace detail
{
void
car2top_matrix(double sinf, double sinl, double cosf, double cosl, double* coef)
noexcept;

void
car2top_cov_matrix(double sin2f, double sin2l, double cos2f, double cos2l,
    double* coef) noexcept;
} // namespace detail

} // namespace ngpt

#endif

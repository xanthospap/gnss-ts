#include <cmath>
#include "geodesy.hpp"

/// Compute the rotation matrix R, to convert a cartesian vector to a
/// topocentric one, i.e.
/// [n,e,u]**T = R*[dx, dy, dz]**T
///
/// \param[in] sinf sin(latitude)
/// \param[in] sinl sin(longtitude)
/// \param[in] cosf cos(latitude)
/// \param[in] cosl cos(longtitude)
/// \param[in] coef At output holds the coefficients of the rotation matrix R,
///                 in row major form; must have size >= 9.
///
void
ngpt::detail::car2top_matrix(double sinf, double sinl, double cosf, double cosl,
    double* coef) noexcept
{
    coef[0] = -sinf*cosl; coef[1] = -sinf*sinl; coef[2] = cosf;
    coef[3] = -sinl;      coef[4] = cosl;       coef[5] = 0.0;
    coef[6] = cosf*cosl;  coef[7] = cosf*sinl;  coef[8] = sinf;
}

/// Compute the rotation matrix R, to convert a cartesian vector of variances
/// to topocentric one, i.e.
/// [sn**2,se**2,su**2]**T = R*[sx**2, sy**2, sz**2]**T
///
/// \param[in] sinf sin**2(latitude)
/// \param[in] sinl sin**2(longtitude)
/// \param[in] cosf cos**2(latitude)
/// \param[in] cosl cos**2(longtitude)
/// \param[in] coef At output holds the coefficients of the rotation matrix R,
///                 in row major form; must have size >= 9.
///
void
ngpt::detail::car2top_cov_matrix(double sin2f, double sin2l, double cos2f,
    double cos2l, double* coef) noexcept
{
    coef[0] = sin2f*cos2l; coef[1] = sin2f*sin2l;  coef[2] = cos2f;
    coef[3] = sin2l;       coef[4] = cos2l;        coef[5] = 0.0;
    coef[6] = cos2f*cos2l; coef[7] = cos2f*sin2l;  coef[8] = sin2f;
}

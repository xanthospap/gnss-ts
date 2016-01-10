#ifndef _CARTESIAN_TO_TOPOCENTRIC_
#define _CARTESIAN_TO_TOPOCENTRIC_

#include "car2ell.hpp"
#include <cmath>

namespace geodesy {

/** \details  Transform a vector in cartesian coordinates to the topocentric,
 *            local system around point i (i.e. North(i), East(i), Up(i)). This
 *            is a template function, depending on the ellipsoid parameter;
 *            see ellipsoid.hpp
 *
 *  \parameter[in]  xi     Cartesian, x-component of point i, meters.
 *  \parameter[in]  yi     Cartesian, y-component of point i, meters.
 *  \parameter[in]  zi     Cartesian, z-component of point i, meters.
 *  \parameter[in]  xj     Cartesian, x-component of point j, meters.
 *  \parameter[in]  yj     Cartesian, y-component of point j, meters.
 *  \parameter[in]  zj     Cartesian, z-component of point j, meters.
 *  \parameter[out] north  Vector north component, meters.
 *  \parameter[out] east   Vector east component, meters.
 *  \parameter[out] up     Vector up component, meters.
 * 
 *  \throw    Does not throw.
 *
 *  \note     The ellispoid is needed to transform the cartesian coordinates of
 *            the (reference) point i to ellispoidal coordinates.
 *
 * Reference: Physical Geodesy, p. 209
 */
template<ELLIPSOID E>
void
car2top(const double& xi,const double& yi, const double& zi,
        const double& xj,const double& yj, const double& zj,
        double& north, double& east, double& up)
noexcept
{

    // Ellipsoidal coordinates of reference point.
    double phi_i,lambda_i,h_i;

    // Cartesian to ellipsoidal for reference point.
    car2ell<E>(xi,yi,zi,phi_i,lambda_i,h_i);

    // Trigonometric numbers.
    double cosf { std::cos(phi_i) };
    double cosl { std::cos(lambda_i) };
    double sinf { std::sin(phi_i) };
    double sinl { std::sin(lambda_i) };

    // Catresian vector.
    double dx { xj - xi };
    double dy { yj - yi };
    double dz { zj - zi };

    // Topocentric vector.
    north = - sinf * cosl * dx - sinf * sinl * dy + cosf * dz;
    east  = - sinl * dx        + cosl * dy;
    up    =   cosf * cosl * dx + cosf * sinl * dy + sinf * dz;

    // Finished.
    return;
}

} // end namespace

#endif

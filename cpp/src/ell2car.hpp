#ifndef __ELLIPSOIDAL_TO_CARTESIAN__
#define __ELLIPSOIDAL_TO_CARTESIAN__

#include <cmath>
#include "ellipsoid.hpp"

namespace ngpt {

/** \details  Transform ellipsoidal, geocentric coordinates to cartesian. This 
 *            is a template function, depending on the ellipsoid parameter; 
 *            see ellipsoid.hpp
 *
 *  \param[in]   phi    Ellipsoidal latitude, radians.
 *  \param[in]   lambda Ellipsoidal longtitude, radians.
 *  \param[in]   h      Ellipsoidal height, meters.
 *  \param[out]  x      Cartesian, x-component, meters.
 *  \param[out]  y      Cartesian, y-component, meters.
 *  \param[out]  z      Cartesian, z-component, meters.
 * 
 *  \throw    Does not throw.
 */
template<ellipsoid E = ellipsoid::wgs84>
void 
ell2car(double phi, double lambda, double h,
        double& x,  double& y,     double& z)
noexcept
{
    // Eccentricity squared.
#if __cplusplus > 201103L
    constexpr
#endif
    double e2 { ngpt::eccentricity_squared<E>() };

    // Trigonometric numbers.
    double sinf { std::sin(phi) };
    double cosf { std::cos(phi) };
    double sinl { std::sin(lambda) };
    double cosl { std::cos(lambda) };

    // Radius of curvature in the prime vertical.
    const double N { ellipsoid_traits<E>::a / 
                     std::sqrt(1.0e0-(e2*sinf)*sinf) };

    // Compute geocentric rectangular coordinates.
    x = (N+h) * cosf * cosl;
    y = (N+h) * cosf * sinl;
    z = ((1.0e0-e2) * N + h) * sinf;

    // Finished.
    return;
}

} // end namespace

#endif

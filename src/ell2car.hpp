#ifndef _ELLIPSOIDAL_TO_CARTESIAN_
#define _ELLIPSOIDAL_TO_CARTESIAN_

#include "ellipsoid.hpp"
#include <cmath>

namespace geodesy {

/** \details  Transform ellipsoidal, geocentric coordinates to cartesian. This 
 *            is a template function, depending on the ellipsoid parameter; 
 *            see ellipsoid.hpp
 *
 *  \parameter[in]   phi    Ellipsoidal latitude, radians.
 *  \parameter[in]   lambda Ellipsoidal longtitude, radians.
 *  \parameter[in]   h      Ellipsoidal height, meters.
 *  \parameter[out]  x      Cartesian, x-component, meters.
 *  \parameter[out]  y      Cartesian, y-component, meters.
 *  \parameter[out]  z      Cartesian, z-component, meters.
 * 
 *  \throw    Does not throw.
 */
template<ELLIPSOID E>
void 
ell2car(const double& phi,const double& lambda,const double& h,
        double& x,double& y,double& z) noexcept
{
    // Eccentricity squared.
    constexpr double e2 { EllipsoidTraits<E>::eccentricitySquared() };

    // Trigonometric numbers.
    double sinf { std::sin(phi) };
    double cosf { std::cos(phi) };
    double sinl { std::sin(lambda) };
    double cosl { std::cos(lambda) };

    // Radius of curvature in the prime vertical.
    const double N { EllipsoidTraits<E>::a / std::sqrt(1.0e0-(e2*sinf)*sinf) };

    // Compute geocentric rectangular coordinates.
    x = (N+h) * cosf * cosl;
    y = (N+h) * cosf * sinl;
    z = ((1.0e0-e2) * N + h) * sinf;

    // Finished.
    return;
}

} // end namespace

#endif

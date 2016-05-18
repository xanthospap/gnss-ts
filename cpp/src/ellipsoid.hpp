#ifndef __REFERENCE_ELLIPSOID__
#define __REFERENCE_ELLIPSOID__

#include <cmath>

namespace ngpt {

/// A list of bulitin ellipsoids.
enum class ellipsoid : char {
    grs80,
    wgs84,
    pz90
};

/// Specify ellipsoid parameters for each enumerated ellipsoid model.
template<ellipsoid E>
    struct ellipsoid_traits {};

/// Specify ellipsoid parameters for ellipsoid::grs80.
template<>
    struct ellipsoid_traits<ellipsoid::grs80>
{
    static constexpr double a      { 6378137.0e0 };
    static constexpr double f      { 1.0e00/298.257222101e0 };
    static constexpr const char* n { "GRS80" };
};

/// Specify ellipsoid parameters for ellipsoid::wgs84.
template<>
    struct ellipsoid_traits<ellipsoid::wgs84>
{
    static constexpr double a      { 6378137.0e0 };
    static constexpr double f      { 1.0e00/298.257223563e0 };
    static constexpr const char* n { "WGS84" };
};

/// Specify ellipsoid parameters for ellipsoid::pz90.
template<>
    struct ellipsoid_traits<ellipsoid::pz90>
{
    static constexpr double a      { 6378135.0e0 };
    static constexpr double f      { 1.0e00/298.257839303e0 };
    static constexpr const char* n { "PZ90" };
};

/// Compute the squared eccentricity for an ellipsoid model (i.e. e**2).
template<ellipsoid E>
    constexpr double eccentricity_squared() noexcept
{
    constexpr double f {ellipsoid_traits<E>::f};
    return (2.0e0-f) * f;
}

/// Compute the semi-minor axis (i.e. b) for an ellipsoid model.
template<ellipsoid E>
    constexpr double semi_minor() noexcept
{
  constexpr double a {ellipsoid_traits<E>::a};
  return a - eccentricity_squared<E>() * a;
}

/// Compute the Normal radious of curvature, for a given latitude.
/// Physical Geodesy, p. 194
template<ellipsoid E>
    double N (double lat) noexcept
{
    constexpr double a {ellipsoid_traits<E>::a};

    double cosf  { std::cos(lat) };
    double sinf  { std::sin(lat) };
    double acosf { a * cosf };
    double bsinf { semi_minor<E>() * sinf };
    double den   { std::sqrt(acosf*acosf + bsinf*bsinf) };

    return (a * a) / den;
}

} // end namespace geodesy

#endif

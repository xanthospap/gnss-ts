#ifndef __REFERENCE_ELLIPSOID__
#define __REFERENCE_ELLIPSOID__

#include <cmath>

namespace geodesy {

enum class ELLIPSOID : char {
    GRS80,
    WGS84,
    PZ90
};

template<ELLIPSOID E>
struct EllipsoidTraits { };

template<>
struct EllipsoidTraits<ELLIPSOID::GRS80>
{
    static constexpr double a      { 6378137.0e0 };
    static constexpr double f      { 1.0e00/298.257222101e0 };
    static constexpr const char* n { "GRS80" };
};

template<>
struct EllipsoidTraits<ELLIPSOID::WGS84>
{
    static constexpr double a      { 6378137.0e0 };
    static constexpr double f      { 1.0e00/298.257223563e0 };
    static constexpr const char* n { "WGS84" };
};

template<>
struct EllipsoidTraits<ELLIPSOID::PZ90>
{
    static constexpr double a      { 6378135.0e0 };
    static constexpr double f      { 1.0e00/298.257839303e0 };
    static constexpr const char* n { "PZ90" };
};

template<ELLIPSOID E>
constexpr double eccentricity_squared()
noexcept
{
  constexpr double f { EllipsoidTraits<E>::f };
  return ( 2.0e0 - f ) * f;
}

template<ELLIPSOID E>
constexpr double semi_minor()
noexcept
{
  constexpr double a { EllipsoidTraits<E>::a };
  return a - eccentricity_squared<E>() * a;
}

/// Normal radious of curvature
/// Physical Geodesy, p. 194
template<ELLIPSOID E>
double N(double lat) noexcept
{
    constexpr double a { EllipsoidTraits<E>::a };

    double cosf  { std::cos(lat) };
    double sinf  { std::sin(lat) };
    double acosf { a * cosf };
    double bsinf { semi_minor<E>() * sinf };
    double den   { std::sqrt(acosf*acosf + bsinf*bsinf) };

    return (a * a) / den;
}

} // end namespace geodesy

#endif

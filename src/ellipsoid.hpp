#ifndef _REFERENCE_ELLIPSOID_
#define _REFERENCE_ELLIPSOID_

# include <cmath>

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
    static constexpr double a { 6378137.0e0 };
    static constexpr double f { 1.0e00/298.257222101e0 };
    static constexpr const char* n { "GRS80" };
    static constexpr double eccentricitySquared() noexcept
    { return ( 2.0e0 - f ) * f; }
    static constexpr double semiMinor() noexcept
    { return a - eccentricitySquared() * a; }
};

template<>
struct EllipsoidTraits<ELLIPSOID::WGS84>
{
    static constexpr double a { 6378137.0e0 };
    static constexpr double f { 1.0e00/298.257223563e0 };
    static constexpr const char* n { "WGS84" };
    static constexpr double eccentricitySquared() noexcept
    { return ( 2.0e0 - f ) * f; }
    static constexpr double semiMinor() noexcept
    { return a - eccentricitySquared() * a; }
};

template<>
struct EllipsoidTraits<ELLIPSOID::PZ90>
{
    static constexpr double a { 6378135.0e0 };
    static constexpr double f { 1.0e00/298.257839303e0 };
    static constexpr const char* n { "PZ90" };
    static constexpr double eccentricitySquared() noexcept
    { return ( 2.0e0 - f ) * f; }
    static constexpr double semiMinor() noexcept
    { return a - eccentricitySquared() * a; }
};

/// Normal radious of curvature
/// Physical Geodesy, p. 194
template<ELLIPSOID E>
double N(double lat) noexcept
{
    double cosf  { std::cos(lat) };
    double sinf  { std::sin(lat) };
    double acosf { E::a * cosf };
    double bsinf { E::semiMinor() * sinf };
    double den   { std::sqrt(acosf*acosf + bsinf*bsinf) };
    return (E::a * E::a) / den;
}


} // end namespace geodesy

#endif

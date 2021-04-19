#ifndef __NGPT_EARTHQUAKE_CATALOGUE2_HPP__
#define __NGPT_EARTHQUAKE_CATALOGUE2_HPP__

///
/// @brief This file defines classes and functions to treat earthquake events
///

#include "ggdatetime/dtcalendar.hpp"
#include "ggeodesy/ellipsoid.hpp"
#include "ggeodesy/geodesy.hpp"
#include "ggeodesy/units.hpp"
#include "ggeodesy/vincenty.hpp"
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace ngpt {

/// @struct earthquake
///
/// @brief A simple class to hold an earthquake event.
struct earthquake {
  ngpt::datetime<ngpt::milliseconds> m_epoch; ///< The datetime it happened
  double m_lon;       ///< Epicenter longtitude (radians)
  double m_lat;       ///< Epicenter latitude (radians)
  double m_depth;     ///< Depth (meters)
  double m_magnitude; ///< The magnitude in (??)

  /// @brief Return the distance of a point on the ellipsoid from the epcenter.
  ///
  /// The distance is computed along the geodesic that connects the two points
  /// on the ellipsoid E. The line is from the epicenter to the given point.
  /// To compute the distance, the function uses the (inverse) Vincenty
  /// algorithm; hence, the forward and backward azimouths are also computed
  /// and returned (as parameters frw_az and bkw_az).
  ///
  /// @tparam     E       The reference ellipsoid (default is
  ///                     ngpt::ellipsoid::wgs84)
  /// @param[in]  lat     The latitide of the point P (radians)
  /// @param[in]  lon     The longtitude of the point P (radians)
  /// @param[out] frw_az  Forward azimouth (i.e. epicenter to P)
  /// @param[out] bkw_az  Backward azimouth (i.e. P to epicenter)
  /// @return             The distance from the epicenter to point P on the
  ///                     ellipsoid.
  template <ellipsoid E = ellipsoid::wgs84>
  double epicenter_distance(double lat, double lon, double &frw_az,
                            double &bkw_az) const {
    return inverse_vincenty<E>(m_lat, m_lon, lat, lon, frw_az, bkw_az, 1e-12);
  }

  /// @brief Concatenate the earthquake elements to  a string.
  ///
  /// The elements (aka instance variables) are joined to a string; the
  /// string follows the convention of the NOA published earthquake catalogue
  /// files, that is:
  /// YYYY OOO DD   HH MM SS.S   LAT     LON     DEPTH(km)  MAGNITUDE
  /// where 'OOO' is the month as 3-char uppercase string, LAT and LON are
  /// given in decimal degrees with a precision of e-2, depth is given in
  /// integer km and magnitude in M with precision 1e-1. Example:
  /// 1964 FEB 24   23 30 25.0   38.90   23.90   10         5.3
  std::string to_string() const noexcept;

}; // earthquake

/// @brief Format a datetime<T> instance based on the NOA catalogue files.
///
/// Given a datetime<T> instance, format it as a string of type:
/// YYYY OOO DD  HH MM SS.S
/// Where 'OOO' is the 3 first chars of the month, in uppercase.
/// @tparam    T  datetime instance resolution
/// @param[in] t  datetime<T> instance to be transformed to string
/// @return       string; the input datetime instance formated as:
///               'YYYY OOO DD HH MM SS.S'
std::string strfdt_as_noa(const ngpt::datetime<ngpt::milliseconds> &t);

} // namespace ngpt

#endif

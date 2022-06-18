#ifndef __NGPT_DATA_POINT_HPP__
#define __NGPT_DATA_POINT_HPP__

#include "genflags.hpp"
#include "ts_flag.hpp"

namespace dso {
/// @class data_point
///
/// A time-series is a series of data points. This data_point class is designed
/// to assist the handling of timeseries. The class itself does very little
/// and is pretty generic. The only limitation is that the F template parameter
/// is an enumeration class that can act as a flag, i.e. dso::flag<F> makes
/// sense and has a default constructor.
/// For example, coordinate time-series, could use: dso::pt_marker (as the F
/// parameter).
/// Each instance of the data_point class, has:
///     - a value
///     - a sigma (i.e. std. deviation) and
///     - a flag (of type dso::flag<F>)
///
/// @see dso::flag template class
///
/// @todo there should be a restriction on F that the function:
///       'bool skip(dso::flag<F>) noexcept' exists.
///
/// @tparam F An enumeration class to act as flag; F must have:
///           - a default constructor
///           - a function with signature 'bool skip(dso::flag<F>) noexcept'
///
struct data_point {
  /// Simplify the flag type.
  using dp_flag = dso::flag<dso::pt_marker>;

  /// Constructor.
  data_point(double val = 0e0, double sigma = 1e-3,
             dp_flag f = dp_flag{}) noexcept
      : m_value{val}, m_sigma{sigma}, m_flag{f} {}

  /// (const) get the value.
  double value() const noexcept { return m_value; }

  /// (const) get the sigma (std. dev).
  double sigma() const noexcept { return m_sigma; }

  /// (const) get the flag.
  dp_flag flag() const noexcept { return m_flag; }

  /// Should the data point be skipped/ignored ?
  /// @todo what the fuck is this? see also the detailed class description.
  /*bool skip() const noexcept { return dso::__skip__(this->m_flag); }*/

  /// equality operator
  bool operator==(const data_point &other) const noexcept {
    return (m_value == other.m_value && m_sigma == other.m_sigma) &&
           (m_flag == other.m_flag);
  }
  bool operator!=(const data_point &other) const noexcept {
    return !(*this == other);
  }

  bool is_clean() const noexcept { return m_flag.is_clean(); }

  double m_value; ///< The data point's value
  double m_sigma; ///< The data point's sigma (i.e. standard deviation)
  dp_flag m_flag; ///< The point's flag

}; // data_point

} // namespace dso
#endif

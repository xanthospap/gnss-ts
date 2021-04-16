#ifndef __TIMESERIES_FLAGS2_HPP__
#define __TIMESERIES_FLAGS2_HPP__

///
/// @file  ts_flag.hpp
///
/// @brief
///

#include <cstdint>

namespace ngpt {

/// @enum pt_marker
///
/// An enumeration type, to hold possible flags for (coordinate) time-series
/// data points.
///
/// @warning If a new pt_marker enum is added (or removed), most of the
///          functions in this file should change!
enum class pt_marker : uint_fast8_t {
  outlier = 0x0001, ///< Signify an outlier
  skip = 0x0002     ///< Signify a data-point that must be skipped
};

struct data_flag {
  uint_fast8_t bitflag{0};

  /// Set a pt_marker.
  /// @param[in] f  A pt_marker instance; this will be "switched on".
  void set(pt_marker f) noexcept { bitflag |= static_cast<uint_fast8_t>(f); }

  /// Clear a pt_marker.
  /// @param[in] f  Clear (i.e. "switch off") this pt_marker.
  void clear(pt_marker f) noexcept {
    bitflag &= ~(static_cast<uint_fast8_t>(f));
  }

  /// Clear the instance from all pt_markers. All pt_marker are "switched off".
  void clear() noexcept { bitflag = static_cast<uint_fast8_t>(0); }

  /// Check if a pt_marker is set (i.e. on).
  /// @param[in] f  Check if this pt_marker is "switched on".
  /// @return true if pt_marker f is set; false otherwise.
  bool is_set(pt_marker f) const noexcept {
    return bitflag & static_cast<uint_fast8_t>(f);
  }

  /// Check if a flag is clean (nothing is set).
  /// @return true if instance has no pt_marker set; false otherwise.
  bool is_clean() const noexcept { return !static_cast<uint_fast8_t>(bitflag); }

  /// Equality operator (two instances have the same pt_marker's on)
  /// @param[in] f Instance to check against
  /// @return true if the two instances have exactly the same pt_marker s set;
  ///         false otherwise.
  bool operator==(data_flag f) const noexcept { return bitflag == f.bitflag; }

  /// InEquality operator.
  /// @param[in] f Instance to check against
  /// @return false if the two instances have exactly the same pt_marker s set;
  ///         else true.
  bool operator!=(data_flag f) const noexcept { return !(this->operator==(f)); }

}; // data_flag

} // namespace ngpt

#endif

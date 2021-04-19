#ifndef __TIMESERIES_FLAGS_HPP__
#define __TIMESERIES_FLAGS_HPP__

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

/// @enum ts_event
///
/// A (strongly typed) enumeration type to hold possible events for coordinate
/// time-series.
///
/// @warning Any changes here (e.g. adding a new ts_event, will affect a big
/// part of the rest of the code (e.g. event, event_list, etc).
enum class ts_event : uint_fast8_t {
  jump = 1,           ///< A jump in the time-series (i.e. offset)
  earthquake = 2,     ///< An earthquake
  velocity_change = 4 ///< A velocity change
};

} // namespace ngpt

#endif

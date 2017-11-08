#ifndef __TIMESERIES_FLAGS_HPP__
#define __TIMESERIES_FLAGS_HPP__

///
/// @file  tsflag.hpp
///
/// @brief This file defines the (strongly-typed) enumeration types, to be used
///        for the timeseries<...> and crdts<...> classes.
///
/// @todo Add an example (somewhere) of how a user can adapt the enums, or make
///       new ones to use for a time-series.

#include <iostream>
#include "genflags.hpp"

namespace ngpt
{

/// An enumeration type, to hold possible flags for (coordinate) time-series 
/// data points.
///
/// @warning If a new pt_marker enum is added (or removed), most of the
///          functions in this file should change!
enum class pt_marker
: int
{
    outlier = 1, ///< Signify an outlier
    skip    = 2  ///< Signify a data-point that must be skipped
};

/// Check if a data-point with a certain flag should be ignored (i.e. skipped).
bool
__skip__(flag<pt_marker> p) noexcept;

/// @brief Write a flag<pt_marker> instance.
///
/// For any enumeration type that can be wrapped around the flag (template)
/// class, there should be an overload for the '<<' operator. I.e. how are
/// we supposed to 'write' this pt_marker(s)?
/// So, this is the convention:
///   - a pt_marker::outlier is written as 'o'
///   - a pt_marker::skip    is written as 's'
/// A clean flag<pt_marker> (i.e. if no pt_marker is set), will write nothing!
/// If a flag<pt_marker> is flagged both as outlier and as skip, then the
/// string 'os' will be written.
///
/// @param[in] os      The stream to write the instance at.
/// @param[in] marker  The flag<pt_marker> to write.
/// @return            The output stream (after writting)
std::ostream&
operator<<(std::ostream& os, const flag<pt_marker>& marker)
{
    if ( marker.check(pt_marker::outlier) ) os << 'o';
    if ( marker.check(pt_marker::skip) ) os << 's';

    return os;
}

/// A (strongly typed) enumeration type to hold possible flags for coordinate
/// time-series events.
enum class ts_event
: int
{
    jump               = 1, ///< A jump in the time-series (i.e. offset)
    earthquake         = 2, ///< An earthquake
    velocity_change    = 4  ///< A velocity change
};

/// @brief Convert a ts_event to its identifing character.
char
event2char(ts_event event) noexcept;

} // namespace ngpt

#endif

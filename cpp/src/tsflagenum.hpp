/// \file tsflagenum.hpp
///
/// \brief Enumeration class to act as flags for time-series.
///
/// \see ngpt::flag
///

#ifndef __TS_FLAG_ENUM__
#define __TS_FLAG_ENUM__

namespace ngpt
{

/// An enumeration type to holf possible flags for coordinate time-series.
enum class ts_events : char
{
    jump,
    earthquake,
    velocity_change,
    outlier,
    skip
};

} // end namespace ngpt

#endif

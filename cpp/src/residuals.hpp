#ifndef __NGPT_TIMESERIES__
#define __NGPT_TIMESERIES__

// standard headers
#include <vector>
#include <algorithm>
#ifdef DEBUG
#include <iostream>
#include <cstdio>
#endif

// Eigen headers
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/QR"

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "timeseries.hpp"

namespace ngpt
{

template<class T,
        class F,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    void
    outlier_detection_nikolaidis(timeseries<T, F>& residuals)
{
}

} // end namespace ngpt

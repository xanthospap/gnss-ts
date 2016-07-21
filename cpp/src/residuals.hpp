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
    class reg_residuals
{
public:
    /// The specific datetime<T> class we will be using.
    using epoch = ngpt::datetime<T>;
    
    /// Simplify the flag type.
    using tflag = ngpt::flag<F>;

    /// The data points
    using entry = ngpt::data_point<F>;
    
    /// An event is described by the event type and a time-stamp (i.e. epoch).
    using event = std::pair<epoch, tflag>;

    explicit reg_residuals(const Eigen::VectorXd& v, const std:vector<epoch>* epoch_ptr,
        double sigma)
    res_ts{v.rows()},
    aposteriori_sigma{sigma}
    {
        std::size_t sz = v.rows();
        for (std::size_t i = 0;i < sz; i++) {

    }

private:
    timeseries<T, F> res_ts;
    double           aposteriori_sigma;
}

} // end namespace ngpt

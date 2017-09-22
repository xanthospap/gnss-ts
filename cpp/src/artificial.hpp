#ifndef __NGPT_ARTFICIAL_TS_HPP__
#define __NGPT_ARTFICIAL_TS_HPP__

// c++ standard headers
#include <stdexcept>
#include <cmath>
#include <random>
#include <iostream>

// gtms headers
#include "crdts.hpp"

namespace ngpt
{

/// \brief Compute harsh number of points, from start to stop with step dt.
///
/// If we want to make a grid of points (i.e. an abscissas) in the interval
/// [start, stop) with step = dt, this function will return a harsh estimate
/// (an overestimate actually) of the number of points this abscissa will
/// contain.
/// This function will actually compute something like:
/// (stop - start) / dt in days.
///
/// \parameter[in] start 
/// \parameter[in] stop
/// \parameter[in] dt
/// \return        Number of days between start and stop, aka
///                (stop - start) / dt in days.
template<class T>
    std::size_t
    pts_in_interval(datetime<T> start, datetime<T> stop, datetime_interval<T> dt)
{
    assert( stop >= start );
    double max_in_day = static_cast<double>( T::max_in_day );
    auto   diff       = stop.delta_date(start);
    double dpts       = static_cast<double>(diff.days())
                        + diff.secs()/max_in_day;
    double ddt        = static_cast<double>(dt.days())
                        + dt.secs()/max_in_day;
    assert( diff>0 && ddt>0 );
    return static_cast<std::size_t>(dpts/ddt+1);
}

/// \brief Construct a synthetic time-series based in input parameters.
///
/// Construct a time-series with its data points following the input model. The
/// time-series will span the period [start, stop) with data-points every dt
/// (i.e. the step = dt).
/// Noise will be added to the data-points, following a normal distribution
/// with mean = mean and std. deviation = stddev.
///
/// \warning This function will allocate memory for the epoch vector; the user
///          is responsible for deallocating this vector!
///
/// \parameter[in] start  (inlusive) start datetime of the time-series.
/// \parameter[in] stop   (non-inlusive) ending datetime of the time-series.
/// \parameter[in] dt     Step size for data-points.
/// \parameter[in] model  A model to characterize the data-points; every data-points
///                       at each instance in the range start:stop:dt , is computed
///                       via this model.
/// \parameter[in] mean   The mean of a Gaussian distribution to be added to
///                       the data-points (default is 0).
/// \parameter[in] stddev The std. deviation of a Gaussian distribution to be 
///                       added to the data-points (default is 0.5).
/// \return               A time-series instance.
///
template<class T, class F>
    timeseries<T,F>
    synthetic_ts(datetime<T> start, datetime<T> stop, datetime_interval<T> dt,
        const ts_model<T>& model, double mean=0e0, double stddev=0.5e0)
{
    std::size_t harsh_size = pts_in_interval(start, stop, dt);

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d(mean, stddev);

    std::vector<datetime<T>>* epochs = new std::vector<datetime<T>>;
    epochs->reserve(harsh_size);
    timeseries<T,F> ts ( epochs );

    for (auto t = start; t < stop; t+=dt) epochs->emplace_back(t);
    auto vec = model.make_model(*epochs);
    for (auto v : vec) ts.add_point( v+d(gen) );

    return ts;
}

/// \brief Construct a synthetic time-series based in input parameters.
///
/// Construct a time-series with its data points following the input model. The
/// time-series will span the period [epochs[0], epochs[1], ... ,epochs[N]].
/// Noise will be added to the data-points, following a normal distribution
/// with mean = mean and std. deviation = stddev.
///
/// \parameter[in] epochs A vector of epochs; for each element, a data-point
///                       value will be computed. This vector, will then be
///                       assigned to the returned timeseries instance.
/// \parameter[in] model  A model to characterize the data-points; every data-points
///                       at each instance in the range start:stop:dt , is computed
///                       via this model.
/// \parameter[in] mean   The mean of a Gaussian distribution to be added to
///                       the data-points (default is 0).
/// \parameter[in] stddev The std. deviation of a Gaussian distribution to be 
///                       added to the data-points (default is 0.5).
/// \return               A time-series instancer; its epochs vector will be
///                       the epochs parameter passed in the function.
///
/// \warning No new epoch vector will be created for the returned timeseries
///          instance; instead, the input parameter 'epochs' will be assigned
///          to the new instance.
///
template<class T, class F>
    timeseries<T,F>
    synthetic_ts(std::vector<datetime<T>>& epochs,
        const ts_model<T>& model, double mean=0e0, double stddev=0.5e0)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d(mean, stddev);

    timeseries<T,F> ts ( &epochs );

    auto vec = model.make_model(epochs);
    for (auto v : vec) ts.add_point( v+d(gen) );

    return ts;
}

}// end ngpt

#endif

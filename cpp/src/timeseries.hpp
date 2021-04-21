#ifndef __NGPT_TIMESERIES_HPP__
#define __NGPT_TIMESERIES_HPP__

#include <algorithm>
#include <ggdatetime/dtfund.hpp>
#include <iterator>
#include <tuple>
#include <vector>
#ifdef DEBUG
#include <cstdio>
#include <iomanip>
#include <iostream>
#endif
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/QR"
#include "ggdatetime/dtcalendar.hpp"
// #include "model.hpp"
#include "genflags.hpp"
#include "ts_flag.hpp"

namespace ngpt {

/// @class data_point
///
/// A time-series is a series of data points. This data_point class is designed
/// to assist the handling of timeseries. The class itself does very little
/// and is pretty generic. The only limitation is that the F template parameter
/// is an enumeration class that can act as a flag, i.e. ngpt::flag<F> makes
/// sense and has a default constructor.
/// For example, coordinate time-series, could use: ngpt::pt_marker (as the F
/// parameter).
/// Each instance of the data_point class, has:
///     - a value
///     - a sigma (i.e. std. deviation) and
///     - a flag (of type ngpt::flag<F>)
///
/// @see ngpt::flag template class
///
/// @todo there should be a restriction on F that the function:
///       'bool skip(ngpt::flag<F>) noexcept' exists.
///
/// @tparam F An enumeration class to act as flag; F must have:
///           - a default constructor
///           - a function with signature 'bool skip(ngpt::flag<F>) noexcept'
///
struct data_point {
  /// Simplify the flag type.
  using dp_flag = ngpt::flag<ngpt::pt_marker>;

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
  /*bool skip() const noexcept { return ngpt::__skip__(this->m_flag); }*/

  /// equality operator
  bool operator==(const data_point &other) const noexcept {
    return (m_value == other.m_value && m_sigma == other.m_sigma) &&
           (m_flag == other.m_flag);
  }
  bool operator!=(const data_point& other) const noexcept {
    return !(*this == other);
  }

  bool is_clean() const noexcept { return m_flag.is_clean(); }

  double m_value; ///< The data point's value
  double m_sigma; ///< The data point's sigma (i.e. standard deviation)
  dp_flag m_flag; ///< The point's flag

}; // data_point

/// @class timeseries
///
/// @brief A generic time-series class.
///
/// @note     A time-series instance does **NOT** own an epoch vector; each
///           instance only holds a pointer to a vector of epochs. The
///           construction/desctruction of this vector, is the responsibility
///           of the user.
///
/// @todo     - Should the time-series be always in correct time-order. Say more
///             about this ....
///           - is the mean values really valid, always ??
///           - using access to data_points (& epochs) via the iterator class,
///             can change the time-serie's mean value etc, without the instance
///             knowing it (that affects mostly the mean value and m_skipped)
///
/// @warning  Mean value and number of skipped points should always be correct
///           (i.e. updated).
///
/// @example  test_ts.cpp
///
class timeseries {
public:
  using dp_flag = ngpt::flag<ngpt::pt_marker>;

  /// Constructor. If a vector of epochs is passed in, then we know we have
  /// our epochs. In this case, reserve space (memory) for the data points.
  ///
  /// @param[in] epochs A pointer to a vector of ngpt::datetime<T> instances.
  ///
  /// @note Even though no data points are added to the time-series, enough
  ///       space is allocated to hold them (only allocated **not**
  ///       initialized!
  explicit timeseries(std::vector<ngpt::datetime<ngpt::milliseconds>> *epochs =
                          nullptr) noexcept
      : m_epochs(epochs) {
    if (m_epochs)
      m_data.reserve(m_epochs->size());
  }

  /// Constructor. Use this constructor if you don't know exactly how many
  /// elements (i.e. data points) the time-series will have, but you do have
  /// a clue.
  ///
  /// @param[in] size_hint  A hint for (or even better the actual) number of
  ///                       data points in the time-series.
  ///
  /// @note The epoch array will be empty after the construction. Users have
  ///       to set it afterwards.
  explicit timeseries(std::size_t size_hint = 100) noexcept
      : m_epochs(nullptr) {
    m_data.reserve(size_hint);
  }

  /// Copy constructor.
  ///
  /// @warning Note that the epoch vector is not (deep) copied; it is only
  ///          set to point to the same epoch vector as the copied-from time-
  ///          series.
  ///
  /// @param[in] ts    The original time-series to be copied.
  /// @param[in] start The index to start copying from; if not given, it is
  ///                  set to 0 (i.e. from start)
  /// @param[in] end   The index of one-past-the-end to stop copying; if not
  ///                  given, it will be set to ts.m_data.size().
  ///
  /// timeseries<...> ts1 { ... };
  /// timeseries<...> ts2{ts1, 10, 100} will copy to ts2 all ts2 values
  /// between indexes [10,...,99]
  timeseries(const timeseries &ts, std::size_t start = 0, std::size_t end = 0);

  /// Move constructor.
  /// @note The resluting time-serie's epoch vector (pointer), will be set to
  /// (point to) the original (i.e. ts).
  timeseries(timeseries &&ts) noexcept
      : m_epochs(ts.m_epochs), m_data{std::move(ts.m_data)} {}

  /// Assignment operator.
  /// @warning Note that the epoch vector is not (deep) copied; it is only
  ///          set to point to the same epoch vector as the copied-from time-
  ///          series.
  timeseries &operator=(const timeseries &ts) noexcept {
    if (this != &ts) {
      m_epochs = ts.m_epochs;
      m_data = ts.m_data;
    }
    return *this;
  }

  /// Move assignment operator.
  timeseries &operator=(timeseries &&ts) noexcept {
    if (this != &ts) {
      m_epochs = ts.m_epochs;
      m_data = std::move(ts.m_data);
    }
    return *this;
  }

  bool operator==(const timeseries& other) const noexcept {
    return m_epochs == other.m_epochs && m_data == other.m_data;
  }
  bool operator!=(const timeseries& other) const noexcept {
    return !this->operator==(other);
  }

  /// Split a time-series; return two new time-series in the interval:
  /// [0-idx) and [idx-end).
  ///
  /// @return A tuple (pair) containing the two new time-series.
  ///
  /// @todo   What the fuck should i do with the epochs of each sub-timeseries??
  auto split(std::size_t idx) const {
    timeseries left(*this, 0, idx);
    timeseries right(*this, idx);
    return std::make_tuple(std::move(left), std::move(right));
  }

  /// Add a data point; returns the new mean value.
  /// @note   The instance's mean value is updated; so is the number of
  ///         skipped data points (if needed).
  /// @return The updated time-series mean value.
  void add_point(double val, double sigma = 1e-3,
                   dp_flag f = dp_flag{}) noexcept {
    m_data.emplace_back(data_point{val, sigma, f});
  }

  std::size_t mark(ngpt::pt_marker type, ngpt::datetime<ngpt::milliseconds> start = ngpt::datetime<ngpt::milliseconds>::min(), ngpt::datetime<ngpt::milliseconds> end=ngpt::datetime<ngpt::milliseconds>::max()) noexcept;

  /// Compute the mean (i.e. central epoch). This version uses the very
  /// first and last epochs to compute the mean, regardless if they are
  /// marked as unused. The mean epoch is obviously half the distance
  /// between the first and last epochs.
  /*
  ngpt::datetime<ngpt::milliseconds> central_epoch() const noexcept {
    auto delta_dt = ngpt::delta_date(last_epoch(), first_epoch());
    auto central_epoch{first_epoch()};
    central_epoch += (delta_dt / 2);
    return central_epoch;
  }*/

#ifdef DEBUG
  std::size_t size() const noexcept { return m_data.size(); }
  const std::vector<data_point>& data_points_vec() const noexcept { return m_data; }
  std::vector<data_point>& data_points_vec() noexcept { return m_data; }
  const std::vector<ngpt::datetime<ngpt::milliseconds>>& epochs_vec() const noexcept { return *m_epochs; }
#endif

private:
  /// A pointer to a vector of datetime<T> instances.
  std::vector<ngpt::datetime<ngpt::milliseconds>> *m_epochs;
  /// The vector of data points.
  std::vector<data_point> m_data;

}; // timeseries

} //ngpt

#endif

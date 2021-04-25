#ifndef __NGPT_TIMESERIES_HPP__
#define __NGPT_TIMESERIES_HPP__

#include <vector>
#ifdef DEBUG
#include <cstdio>
#include <iomanip>
#include <iostream>
#endif
#include "data_point.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include "model.hpp"

namespace ngpt {

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

  timeseries(const ngpt::ts_model &m,
             std::vector<ngpt::datetime<ngpt::milliseconds>> *epochs,
             double white_noise_std_dev = -1e0) noexcept;

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
  timeseries(const timeseries &ts) noexcept
      : m_epochs(ts.m_epochs), m_data(ts.m_data){};

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

  bool operator==(const timeseries &other) const noexcept {
    return m_epochs == other.m_epochs && m_data == other.m_data;
  }
  bool operator!=(const timeseries &other) const noexcept {
    return !this->operator==(other);
  }

  /// Split a time-series; return two new time-series in the interval:
  /// [0-idx) and [idx-end).
  ///
  /// @return A tuple (pair) containing the two new time-series.
  ///
  /// @todo   What the fuck should i do with the epochs of each sub-timeseries??
  /*
  auto split(std::size_t idx) const {
    timeseries left(*this, 0, idx);
    timeseries right(*this, idx);
    return std::make_tuple(std::move(left), std::move(right));
  }
  */

  /// Add a data point; returns the new mean value.
  /// @note   The instance's mean value is updated; so is the number of
  ///         skipped data points (if needed).
  /// @return The updated time-series mean value.
  void add_point(double val, double sigma = 1e-3,
                 dp_flag f = dp_flag{}) noexcept {
    m_data.emplace_back(data_point{val, sigma, f});
  }

  /// Mark all data_points in the m_data vector as type, for the range
  /// [start, end).
  /// parameter[in] start Start of data range (inclusive)
  /// parameter[in] end   Stop marking at this date (non-inclusive)
  /// parameter[in] type  ngpt::pt_marker to set for the data points in the
  ///                     range
  /// @return Number of data points that were marked
  std::size_t mark(ngpt::pt_marker type,
                   ngpt::datetime<ngpt::milliseconds> start =
                       ngpt::datetime<ngpt::milliseconds>::min(),
                   ngpt::datetime<ngpt::milliseconds> end =
                       ngpt::datetime<ngpt::milliseconds>::max()) noexcept;

  /// Compute the mean (i.e. central epoch). This version uses the very
  /// first and last epochs to compute the mean, regardless if they are
  /// marked as unused. The mean epoch is obviously half the distance
  /// between the first and last epochs.
  ngpt::datetime<ngpt::milliseconds> central_epoch() const noexcept;

  /// Cut (in-place) the m_data vector to fit the interval [start, end). The
  /// m_epochs vector WILL NOT BE AFFECTED!. If you want to also cut the epochs
  /// vector, make a copy and pass it in via the epoch_vec parameter.
  /// Example usage:
  /// timeseries ts1( ... ); // some timeseries instance
  /// auto t1 = ngpt::datetime( ... );
  /// auto t2 = ngpt::datetime( ... );
  /// auto ts2 = ts1; // deep copy m_data and shallow copy m_epochs
  /// vector<data_point> copy = ts1.epochs_vec(); // copy the ts1 m_epochs
  /// vector auto sz = ts2.cut(t1, t2, &copy); ts2.set_epoch_vec(&copy); //
  /// don't forget to set a valid epochs vector
  /// @see test/test_timeseries_1.cc
  ///
  /// @param[in] start First date to include in the new timeseries
  /// @param[in] end   Date limit; no date >= to end will be included (aka this
  ///                  is a non-inclusive limit).
  /// @param[in] epoch_vec A (deep) copy of the instance's m_epochs vector;
  ///                  this vector will be cut to only hold dates included in
  ///                  the [start, end) range.
  /// @return The size of the new m_data vector; if something goes wrong a
  ///         negative integer is returned.
  long cut(ngpt::datetime<ngpt::milliseconds> start,
           ngpt::datetime<ngpt::milliseconds> end,
           std::vector<ngpt::datetime<ngpt::milliseconds>> *epoch_vec =
               nullptr) noexcept;

  void
  set_epoch_vec(std::vector<ngpt::datetime<ngpt::milliseconds>> *e) noexcept {
    m_epochs = e;
  }

#ifdef DEBUG
  std::size_t size() const noexcept { return m_data.size(); }
  const std::vector<data_point> &data_points_vec() const noexcept {
    return m_data;
  }
  std::vector<data_point> &data_points_vec() noexcept { return m_data; }
  const std::vector<ngpt::datetime<ngpt::milliseconds>> &
  epochs_vec() const noexcept {
    return *m_epochs;
  }
#endif

private:
  /// A pointer to a vector of datetime<T> instances.
  std::vector<ngpt::datetime<ngpt::milliseconds>> *m_epochs;
  /// The vector of data points.
  std::vector<data_point> m_data;

}; // timeseries

} // namespace ngpt

#endif

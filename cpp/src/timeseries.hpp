#ifndef __NGPT_TIMESERIES__
#define __NGPT_TIMESERIES__

#include <vector>
#include "dtcalendar.hpp"
#include "genflags.hpp"
#include "tsflagenum.hpp"

namespace ngpt
{

/// A time-series is a series of data points.
class data_point {
public:
    /// Simplify the flag type.
    using tflag = ngpt::flag<ngpt::ts_events>;
    
    /// Constructor.
    explicit data_point(double val=0.0, double sigma=1.0, tflag f=tflag{})
    noexcept
    : m_value{val}, m_sigma{sigma}, m_flag{f}
    {}

    ///
    double value() const noexcept { return m_value; }
    double& value() noexcept { return m_value; }

    ///
    double sigma() const noexcept { return m_sigma; }
    double& sigma() noexcept { return m_sigma; }

    ///
    tflag flag() const noexcept { return m_flag; }
    tflag& flag() noexcept { return m_flag; }

private:
    double m_value;
    double m_sigma;
    tflag  m_flag;

}; // end class data_point

/// A generic time-series class
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class timeseries
{
public:
    /// The specific datetime<T> class we will be using.
    using epoch = ngpt::datetime<T>;
    
    /// Simplify the flag type.
    using tflag = ngpt::flag<ngpt::ts_events>;

    /// Constructor.
    explicit timeseries(std::vector<epoch>* epochs=nullptr) noexcept
    : m_epochs(epochs), m_mean{0.0}, m_events{0}, m_skiped{0}
    {
        if ( m_epochs ) {
            m_data.reserve(m_epochs->size());
        }
    }

    /// Get the (pointer to) epoch vector.
    std::vector<epoch>*& epochs() noexcept { return m_epochs; }

    /// Get the (pointer to) epoch vector.
    const std::vector<epoch>* epochs() const noexcept { return m_epochs; }

    /// Get the data point at index i.
    data_point operator[](std::size_t i) const { return m_data[i]; }

    /// Get the data point at index i.
    data_point& operator[](std::size_t i) { return m_data[i]; }

    /// Get the mean value
    double mean() const noexcept { return m_mean; }

    /// Get the size
    std::size_t size() const noexcept { return m_data.size(); }

    /// Get the number of parameters.
    std::size_t events() const noexcept { return m_events; }

    /// Copy constructor. Note that the epoch vector is set to nullptr.
    timeseries(const timeseries& ts, std::size_t start=0, std::size_t end=0)
    : m_epochs(nullptr), m_mean{ts.m_mean}, m_events{ts.m_events}, m_skiped{ts.m_skiped}
    {
        if (start || end) {
            if (start && !end) {
                end = ts.m_data.size();
            }
            if (end < start || end >= ts.m_data.size()) {
                throw std::domain_error("Invalid start/stop indexes for ts split()");
            }
            m_data.reserve(end-start);
            double sz;
            m_mean = 0.0;
            for (std::size_t i = start; i < end; ++i) {
                sz = i - start;
                m_data.emplace_back(ts[i]);
                m_mean = (ts[i].value() + sz*m_mean)/(sz+1.0);
            }
        } else {
            m_data = ts.m_data;
        }
    }

    /// Move constructor. Note that the epoch vector is set to nullptr.
    timeseries(timeseries&& ts) noexcept
    : m_epochs(nullptr), m_mean{std::move(ts.m_mean)},
      m_data{std::move(ts.m_data)},
      m_events{std::move(ts.m_events)},
      m_skiped{std::move(ts.m_skiped)}
    {}

    /// Assignment operator. Note that the epoch vector is set to nullptr.
    timeseries& operator=(const timeseries& ts) noexcept
    {
        if (this != &ts) {
            m_epochs = nullptr;
            m_mean   = ts.m_mean;
            m_data   = ts.m_data;
            m_events = ts.m_events;
            m_skiped = ts.m_skiped;
        }
        return *this;
    }
    
    /// Move assignment operator. Note that the epoch vector is set to nullptr.
    timeseries& operator=(timeseries&& ts) noexcept
    {
        if (this != &ts) {
            m_epochs = nullptr;
            m_mean   = std::move(ts.m_mean);
            m_data   = std::move(ts.m_data);
            m_events = std::move(ts.m_events);
            m_skiped = std::move(ts.m_skiped);
        }
        return *this;
    }

    /// Split a time-series; return two new time-series in the interval:
    /// [0-idx) and [idx-end). Note that the epoch vector s left as is.
    auto split(std::size_t idx) const
    {
        timeseries left  (*this, 0, idx);
        timeseries right (*this, idx);
        return std::make_tuple(std::move(left), std::move(right));
    }

    /// Add a data point; returns the new mean value.
    double add_point(double val, double sigma=1.0, tflag f=tflag{})
    {
        double sz = static_cast<double>(m_data.size());
        m_data.emplace_back(val, sigma, f);
        m_mean = (val + sz*m_mean)/(sz+1.0);
        if (   f.check(ts_events::jump)
            || f.check(ts_events::velocity_change)
            || f.check(ts_events::earthquake) )
        {
            ++m_events;
        }
        if ( f.check(ts_events::outlier)
            || f.check(ts_events::skip) )
        {
            ++m_skiped;
        }
        return m_mean;
    }

    /// Mark a data point given its index.
    void mark(std::size_t index, ts_events f)
    {
        m_data[index].flag().set(f);
        if (   f == ts_events::jump
            || f == ts_events::velocity_change
            || f == ts_events::earthquake )
        {
            ++m_events;
        }
        if ( f == ts_events::outlier 
            || f == ts_events::skip )
        {
            ++m_skiped;
        }
    }

    /// Solve the least squares via QR (@Eigen)
    ???
    qr_ls_solve(std::vector<double>* phases = nullptr)
    {
        if ( !m_epochs ) {
            throw 1

        /// number of cols/parameters = events + a0 + b0 + 2*(periodic_terms)
        std::size_t parameters = events() + 1  + 1  + (phases ? 2*phases->size() : 0);
        /// number of rows/observations = size - (outliers + skiped)
        std::size_t observations = m_data.size() - m_skiped;
        /// indexes
        std::size_t idx = 0;

        for (const auto& it = m_data.cbegin(); it!= m_data.cend(); ++it)
        {

        }
    }

private:
    /// A pointer to a vector of datetime<T> instances.
    std::vector<epoch>* m_epochs;
    /// The average of the data points.
    double m_mean;
    /// The vector of data points.
    std::vector<data_point> m_data;
    /// number of parameters (offsets + vel. changes + earthquakes)
    std::size_t m_events;
    /// number of outliers/skipped points.
    std::size_t m_skiped;

}; // class timeseries

} // end namespace ngpt
#endif

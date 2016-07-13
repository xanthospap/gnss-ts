#ifndef __NGPT_TIMESERIES__
#define __NGPT_TIMESERIES__

// standard headers
#include <vector>

// Eigen headers
#ifdef KOKO
    #include "Eigen/Core"
    #include "Eigen/QR"
#endif

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

#include "genflags.hpp"
#include "tsflagenum.hpp"

namespace ngpt
{

/// A time-series is a series of data points. This data_point class is designed
/// to assist the handling of timeseries. The class itself does very little 
/// and is pretty generic. The only limitation is that the F template parameter
/// is a enumeration class that can act as a flag, i.e. ngpt::flag<F> makes sense
/// and has a default constructor.
/// For example, coordinate time-series, could use: ngpt::ts_events (as the F
/// parameter).
///
/// \see ngpt::flag template class
///
template<class F> class data_point
{
public:

    /// Simplify the flag type.
    using tflag = ngpt::flag<F>;
    
    /// Constructor.
    explicit data_point(double val=0.0, double sigma=1.0, tflag f=tflag{})
    noexcept
    : m_value{val},
      m_sigma{sigma},
      m_flag{f}
    {}

    /// const get.
    double value() const noexcept { return m_value; }

    /// get/set
    double& value() noexcept { return m_value; }

    /// const get
    double sigma() const noexcept { return m_sigma; }

    /// get/set
    double& sigma() noexcept { return m_sigma; }

    /// const get
    tflag flag() const noexcept { return m_flag; }

    /// get/set
    tflag& flag() noexcept { return m_flag; }

private:
    double m_value; ///< The data point's value
    double m_sigma; ///< The data point's sigma (i.e. standard deviation)
    tflag  m_flag;  ///< The point's flag

}; // end class data_point

/// A generic time-series class
template<class T,
        class F,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class timeseries
{
public:
    /// The specific datetime<T> class we will be using.
    using epoch = ngpt::datetime<T>;
    
    /// Simplify the flag type.
    using tflag = ngpt::flag<F>;

    /// Constructor. If a vector of epochs is passed in, then we know we have
    /// our epochs.
    explicit timeseries(std::vector<epoch>* epochs=nullptr) noexcept
    : m_epochs(epochs),
      m_mean{0.0},
      m_skiped{0}
    {
        if ( m_epochs ) m_data.reserve(m_epochs->size());
    }
    
    /// Constructor. Use this constructor if we don't know how many elements the
    /// time-series will have, but we do have a clue.
    explicit timeseries(std::size_t size_hint) noexcept
    : m_epochs(nullptr),
      m_mean{0.0}
      m_skiped{0}
    {
        m_data.reserve(size_hint);
    }

    /// Get the (pointer to) epoch vector.
    std::vector<epoch>*& epochs() noexcept { return m_epochs; }

    /// Get the (pointer to) epoch vector (const version).
    const std::vector<epoch>* epochs() const noexcept { return m_epochs; }

    /// Get the data point at index i (const version).
    data_point operator[](std::size_t i) const { return m_data[i]; }

    /// Get the data point at index i.
    data_point& operator[](std::size_t i) { return m_data[i]; }

    /// Get the mean value
    double mean() const noexcept { return m_mean; }

    /// Get the size
    std::size_t size() const noexcept { return m_data.size(); }

    /// Get the number of parameters.
    std::size_t events() const noexcept { return m_events; }

    /// Get the first epoch
    /// FIXME what if the first data point is flaged outlier or skipped ??
    epoch first_epoch() const noexcept { return m_epochs[0]; }

    /// Get the last epoch
    epoch last_epoch() const noexcept { return m_epochs[m_epochs.size()-1]; }

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
        /// FIXME if this point is already marked do NOT augment the counter
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

    /// Compute the mean (i.e. central epoch).
    epoch
    central_epoch() const noexcept
    {
        auto delta_dt = ngpt::delta_date(first_epoch(), last_epoch());
        auto central_dt = m_epochs[0].add(std::get<0>(delta_dt),
            std::get<1>(delta_dt));
        return central_epoch;
    }

#ifdef KOKO
    /// Solve the least squares via QR (@Eigen)
    ???
    qr_ls_solve(std::vector<double>* phases = nullptr)
    {
        if ( !m_epochs ) { throw 1 }

        /// number of cols/parameters = events + a0 + b0 + 2*(periodic_terms)
        std::size_t parameters = events() + 1  + 1  + (phases ? 2*phases->size() : 0);
        /// number of rows/observations = size - (outliers + skiped)
        std::size_t observations = m_data.size() - m_skiped;
        /// indexes
        std::size_t idx{0}, counter{0};

        if ( !parameters ) { throw 1; }
        if ( observations < parameters ) { throw 1; }
        
        Eigen::MatrixXd A = Eigen::MatrixXd(observations, parameters);
        Eigen::VectorXd b = Eigen::VectorXd(observations);

        // TODO
        // This is cool for a day interval but probably not enough for more
        // dense sampling rates.
        double mean_epoch { this->central_epoch().as_mjd() };
        for (const auto& it = m_data.cbegin(); it!= m_data.cend(); ++it)
        {
            if ( !it->flag().check(ts_event::outlier)
                && !it->flag().check(ts_event::skip) )
            {
                A.coeff(idx, 0) = 1.0e0;
                A.coeff(idx, 1) = m_epochs->[counter] - mean_epoch;
                ++idx;
            }
            ++counter;
        }
        ///FIXME
    }
#endif

private:
    /// A pointer to a vector of datetime<T> instances.
    std::vector<epoch>* m_epochs;
    /// The average value of the data points.
    double m_mean;
    /// The vector of data points.
    std::vector<data_point> m_data;
    /// Number of outliers/skipped points.
    std::size_t m_skiped;

}; // class timeseries

} // end namespace ngpt
#endif

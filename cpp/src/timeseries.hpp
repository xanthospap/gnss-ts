#ifndef __NGPT_TIMESERIES__
#define __NGPT_TIMESERIES__

// standard headers
#include <vector>
#include <algorithm>

// Eigen headers
#ifdef KOKO
    #include "Eigen/Core"
    #include "Eigen/QR"
#endif

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "genflags.hpp"

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
/// \todo there should be a restriction on F that the function: bool skip(ngpt::flag<F>) noexcept
///       exists.
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

    /// Should the data point be skipped/ignored?
    bool skip() const noexcept { return /*ngpt::__skip__(this->m_flag);*/false; }

private:
    double m_value; ///< The data point's value
    double m_sigma; ///< The data point's sigma (i.e. standard deviation)
    tflag  m_flag;  ///< The point's flag

}; // end class data_point

/// A generic time-series class
/// Mean value and number of skipped points should always be correct (i.e. updated).
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

    /// The data points
    using entry = ngpt::data_point<F>;

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
      m_mean{0.0},
      m_skiped{0}
    {
        m_data.reserve(size_hint);
    }

    /// Get the (pointer to) epoch vector.
    std::vector<epoch>*& epochs() noexcept { return m_epochs; }

    /// Get the (pointer to) epoch vector (const version).
    const std::vector<epoch>* epochs() const noexcept { return m_epochs; }

    /// Get the data point at index i (const version).
    entry operator[](std::size_t i) const { return m_data[i]; }

    /// Get the data point at index i.
    entry& operator[](std::size_t i) { return m_data[i]; }

    /// Get the mean value
    double mean() const noexcept { return m_mean; }

    /// Get the number of data points.
    std::size_t size() const noexcept { return m_data.size(); }

    /// Get the first epoch
    epoch first_epoch() const noexcept { return m_epochs[0]; }

    /// Get the first epoch **NOT** skipped. Returns the value of the first,
    /// valid epoch and sets the idx parameter to its index.
    ///
    /// \todo Should i allow this to throw?
    epoch first_valid_epoch(std::size_t idx) const noexcept
    {
        auto it = std::find_if(std::cbegin(*m_epochs), std::cend(*m_epochs),
            [](const entry& i){return !(i.skip());});
        idx = std::distance(std::cbegin(*m_epochs), it);
        return *it;
    }

    /// Get the last epoch **NOT** skipped. Returns the value of the last,
    /// valid epoch and sets the idx parameter to its index.
    ///
    /// \todo Should i allow this to throw?
    epoch last_valid_epoch(std::size_t idx) const noexcept
    {
        auto it = std::find_if(std::crbegin(*m_epochs), std::crend(*m_epochs),
            [](const entry& i){return !(i.skip());});
        idx = std::distance(it, std::crend(*m_epochs)) - 1;
        return *it;
    }

    /// Get the last epoch
    epoch last_epoch() const noexcept { return m_epochs[m_epochs.size()-1]; }

    /// Copy constructor. Note that the epoch vector is set to nullptr.
    timeseries(const timeseries& ts, std::size_t start=0, std::size_t end=0)
    : m_epochs(nullptr),
      m_mean{ts.m_mean},
      m_skiped{ts.m_skiped}
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
            m_mean   = 0.0;
            m_skiped = 0;
            for (std::size_t i = start; i < end; ++i) {
                sz = i - start;
                m_data.emplace_back(ts[i]);
                m_mean = (ts[i].value() + sz*m_mean)/(sz+1.0);
                if ( m_data[i].skip() ) ++m_skiped;
            }
        } else {
            m_data = ts.m_data;
        }
    }

    /// Move constructor. Note that the epoch vector is set to nullptr.
    timeseries(timeseries&& ts) noexcept
    : m_epochs(nullptr),
      m_mean{std::move(ts.m_mean)},
      m_data{std::move(ts.m_data)},
      m_skiped{std::move(ts.m_skiped)}
    {}

    /// Assignment operator. Note that the epoch vector is set to nullptr.
    timeseries& operator=(const timeseries& ts) noexcept
    {
        if (this != &ts) {
            m_epochs = nullptr;
            m_mean   = ts.m_mean;
            m_data   = ts.m_data;
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
            m_skiped = std::move(ts.m_skiped);
        }
        return *this;
    }

    /// Split a time-series; return two new time-series in the interval:
    /// [0-idx) and [idx-end). Note that the epoch vector is left as is.
    auto split(std::size_t idx) const
    {
        timeseries left  (*this, 0, idx);
        timeseries right (*this, idx);
        return std::make_tuple(std::move(left), std::move(right));
    }


    /// Add a data point; returns the new mean value.
    double add_point(entry&& e)
    {
        double sz = static_cast<double>(m_data.size());
        m_data.emplace_back(e);
        m_mean = (e.value() + sz * m_mean)/(sz+1.0);
        if ( e.skip() ) ++m_skiped;
        return m_mean;
    }

    /// Add a data point; returns the new mean value.
    double add_point(double val, double sigma=1.0, tflag f=tflag{})
    {
        return this->add_point( entry{val, sigma, f} );
    }

    /// Mark a data point given its index.
    void mark(std::size_t index, ts_event f)
    {
        tflag previous = m_data[index].flag();
        m_data[index].flag().set(f);
        if ( !skip(previous) && skip(tflag{f}) ) ++m_skiped;
    }

    /// Compute the mean (i.e. central epoch)
    /*
    epoch
    central_epoch() const noexcept
    {
        auto delta_dt = ngpt::delta_date(first_epoch(), last_epoch());
        auto central_dt = m_epochs[0].add(std::get<0>(delta_dt),
            std::get<1>(delta_dt));
        return central_epoch;
    }
    */

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
    std::vector<entry> m_data;
    /// Number of outliers/skipped points.
    std::size_t m_skiped;

}; // class timeseries

} // end namespace ngpt
#endif

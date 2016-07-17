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
    std::vector<epoch>*& epoch_ptr() noexcept { return m_epochs; }

    /// Get the (pointer to) epoch vector (const version).
    const std::vector<epoch>* epoch_ptr() const noexcept { return m_epochs; }

    /// Get the data point at index i (const version).
    entry operator[](std::size_t i) const { return m_data[i]; }

    /// Get the data point at index i.
    entry& operator[](std::size_t i) { return m_data[i]; }

    /// Get the mean value
    double mean() const noexcept { return m_mean; }

    /// Get the number of data points.
    std::size_t size() const noexcept { return m_data.size(); }

    /// Get the number of epochs.
    std::size_t epochs() const noexcept { return m_epochs ? m_epochs->size() : 0; }

    /// Get the first epoch
    epoch first_epoch() const noexcept { return (*m_epochs)[0]; }
    
    /// Get the last epoch
    epoch last_epoch() const noexcept { return (*m_epochs)[m_epochs->size()-1]; }

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

    /// \brief Compute the mean (i.e. central epoch)
    ///
    /// This version uses the very first and last epochs to compute the mean,
    /// regardless if they are marked as unused.
    ///
    /// \see first_epoch
    /// \see last_epoch
    /// \see cenral_valid_epoch
    epoch
    central_epoch() const noexcept
    {
        using K = typename epoch::sec_type;
        auto delta_dt   = ngpt::delta_date(last_epoch(), first_epoch());
        auto half_lag   = std::get<0>(delta_dt).as_underlying_type() / 2;
        auto sec_of_day = std::get<1>(delta_dt) + 
            ( (std::get<0>(delta_dt).as_underlying_type()%2) ? K{K::max_in_day/2} : K{0} );
        auto central_epoch = (*m_epochs)[0].add(modified_julian_day{half_lag},
            sec_of_day);
        return central_epoch;
    }
    
    /// \brief Compute the mean (i.e. central epoch)
    ///
    /// This version uses the first and last epochs that are not marked as
    /// unused to compute the mean epoch. I.e., if the first 10 epochs --data
    /// points-- are marked as outliers, then they shall not be used to compute
    /// the mean.
    /// 
    /// \see first_valid_epoch
    /// \see last_valid_epoch
    /// \see cenral_epoch
    epoch
    central_valid_epoch() const noexcept
    {
        using K = typename epoch::sec_type;
        auto start_epoch = this->first_valid_epoch();
        auto delta_dt   = ngpt::delta_date(last_valid_epoch(), start_epoch);
        auto half_lag   = std::get<0>(delta_dt).as_underlying_type() / 2;
        auto sec_of_day = std::get<1>(delta_dt) + 
            ( (std::get<0>(delta_dt).as_underlying_type()%2) ? K{K::max_in_day/2} : K{0} );
        auto central_epoch = start_epoch.add(modified_julian_day{half_lag},
            sec_of_day);
        return central_epoch;
    }

    /// Solve the least squares via QR (@Eigen)
    /// To set e.g. a period of 1 year, set periods[0] = 365.25
    /// For 6-months period, periods[0] = 365.25/2
    // Reformulate the problem : A * x = b, with W**2 = P to
    //                          (W*A) * x = (W*b), where W*A=N and W*b=y
    auto
    qr_ls_solve(std::vector<double>* periods = nullptr, double sigma0 = .001)
    {
        if ( !m_epochs ) { throw 1; }
        assert( epochs() == size() );

        /// number of cols/parameters = events + a0 + b0 + 2*(periodic_terms)
        std::size_t parameters = /*m_events.size() +*/ 1  + 1  + (periods ? 2*periods->size() : 0);

        /// number of rows/observations = size - (outliers + skiped)
        std::size_t observations = m_data.size() - m_skiped;
        
        /// indexes
        std::size_t idx{0}, counter{0}, col{0};

        /// set the phases right (trnaform to omegas: 2 * pi * frequency)
        std::vector<double> omegas;
        double freq;
        if ( periods ) {
            omegas = *periods;
            for (auto it = omegas.begin(); it != omegas.end(); ++it) {
                freq = 1.0e0 / *it;
                *it = D2PI * freq;
            }
        }

        if ( !parameters ) { throw 1; }
        if ( observations < parameters ) { throw 2; }
        
        Eigen::MatrixXd A = Eigen::MatrixXd(observations, parameters);
        Eigen::VectorXd b = Eigen::VectorXd(observations);
        Eigen::VectorXd x = Eigen::VectorXd(parameters);

        // TODO
        // This is cool for a day interval but probably not enough for more
        // dense sampling rates.
        double mean_epoch { this->central_epoch().as_mjd() },
               dt, weight;
        // \todo would it be better to fill this column-wise??
        // Instead of forming A and b, we will form A*sqrt(P) and b*sqrt(P)
        for (auto it = m_data.cbegin(); it!= m_data.cend(); ++it)
        {
            if ( !it->skip() ) {
                // delta days from central epoch
                dt = (*m_epochs)[counter].as_mjd() - mean_epoch;
                // weight of observation
                weight = sigma0 / m_data[counter].sigma();

                // coef for constant (linear) term
                A(idx, col) = 1.0e0 * weight;
                ++col;
                // coef for constant (linear) velocity i.e. m/year
                A(idx, col) = weight * (dt / 365.25);
                ++col;
                // Harmonic coefficients for each period ...
                for (auto j = omegas.cbegin(); j != omegas.cend(); ++j) {
                    // cosinus or phase
                    A(idx, col) = std::cos((*j) * dt) * weight;
                    ++col;
                    // sinus or out-of-phase
                    A(idx, col) = std::sin((*j) * dt) * weight;
                    ++col;
                }
                // observation matrix (vector)
                b(idx) = m_data[counter].value() * weight;
                ++idx;
                col = 0;
            }
            ++counter;
        }
        std::cout<<"\nNumber of data points: " << size();
        std::cout<<"\nNumber of epochs:      " << epochs();
        std::cout<<"\nA : "<<observations<<" * "<<parameters;
        std::cout<<"\nb : "<<observations;
        std::cout<<"\nStart: "<<first_epoch().stringify();
        std::cout<<" Stop: "<<last_epoch().stringify();
        std::cout<<" Mean: "<<central_epoch().stringify();
        
        // Solve via QR
        x = A.colPivHouseholderQr().solve(b);
        std::cout<<"\nLS solution vector: \n";
        for (std::size_t i=0;i<parameters;i++) printf("%+15.5f\n",x(i));

        // residual vector u = A*x - b
        Eigen::VectorXd u = Eigen::VectorXd(observations);
        u = A * x - b;

        return x;
    }

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

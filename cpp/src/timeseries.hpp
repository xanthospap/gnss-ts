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

    /// Should the data point be skipped/ignored ?
    bool skip() const noexcept { return ngpt::__skip__(this->m_flag); }

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
    
    /// An event is described by the event type and a time-stamp (i.e. epoch).
    using event = std::pair<epoch, tflag>;

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
        if ( e.skip() ) {
            ++m_skiped;
        } else {
            m_mean = (e.value() + sz * m_mean)/(sz+1.0);
        }
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

    void
    make_model_line(
        const epoch& from,
        const epoch& to,
        const epoch& mean_epoch,
        const std::vector<double>& parameters,
        const std::vector<epoch>& jumps,
        const std::vector<epoch>& vel_changes,
        const std::vector<double>& periods,
        std::vector<double>& model_x,
        std::vector<double>& model_y

    )
    {
        double mindt = from.as_mjd() - mean_epoch.as_mjd(); // in days
        double maxdt = to.as_mjd() - mean_epoch.as_mjd();   // in days
        
        // set the phases right (transform to angular frequencies i.e. omegas:
        //  2 * pi * frequency)
        std::vector<double> omegas {periods};
        double freq;
        for (auto it = omegas.begin(); it != omegas.end(); ++it) {
            freq = 1.0e0 / *it;
            *it = D2PI * freq;
        }

        std::vector<double> x, y;
        x.reserve(maxdt - mindt + 2);
        y.reserve(maxdt - mindt + 2);
        epoch current_epoch{from};
        std::size_t col{0};

        for (double t = mindt; t <= maxdt; t += 1) {
            x.push_back(current_epoch.as_mjd());
            double fx = parameters[0] + parameters[1]*(t/365.25);
            col = 2;
            // Harmonic coefficients for each period ...
            for (auto j = omegas.cbegin(); j != omegas.cend(); ++j) {
                fx += ( parameters[col] * std::cos((*j) * t) +
                    parameters[col+1] * std::sin((*j) * t) );
                col += 2;
            }
            // Set up jumps ...
            for (auto j = jumps.cbegin(); j != jumps.cend(); ++j) {
                if ( *j >= current_epoch ) {
                    fx += parameters[col];
                }
                ++col;
            }
            // Set up velocity changes ...
            for (auto j = vel_changes.cbegin(); j != vel_changes.cend(); ++j) {
                if ( *j >= current_epoch ) {
                    fx += parameters[col]*(t/365.25);
                }
                ++col;
            }
            y.push_back(fx);
            current_epoch = current_epoch.add(modified_julian_day{1}, T{0});
        }
        model_x = std::move(x);
        model_y = std::move(y);
        return;
    }

    // change paramters from pointers to refs.
    auto
    qr_ls_solve(
        const std::vector<epoch>&  jumps,
        const std::vector<epoch>&  vel_changes,
        const std::vector<double>& periods,
        double sigma0 = 1e-03
        /*,std::string* print_model_to=nullptr*/
    )
    {
        assert( m_epochs != nullptr && epochs() == size() );

        // Number of parameters to be estimated:
        std::size_t parameters = 1 + 1              // linear velocity terms
                               + jumps.size()       // jumps/offsets
                               + vel_changes.size() // velocity changes
                               + periods.size()*2;  // harmonic terms

        // Number of observations (ommiting the ones to be skipped)
        std::size_t observations = m_data.size() - m_skiped;
        
        // Can we estimate ?
        if ( observations < parameters ) {
            throw std::runtime_error("[ERROR] Too few observations to perform LS"
            "\n\tNumber of parameters:   "+std::to_string(parameters)+"\n\tNumber of observations: "+std::to_string(observations));
        }

        // set the phases right (transform to angular frequencies i.e. omegas:
        //  2 * pi * frequency)
        std::vector<double> omegas;
        double freq;
        if ( !periods.empty() ) {
            omegas = periods;
            for (auto it = omegas.begin(); it != omegas.end(); ++it) {
                freq = 1.0e0 / *it;
                *it = D2PI * freq;
            }
        }
        
        // Set up the matrices; the model is: A * x = b
        // since we have uncorrelated variance matric though, we will actually
        // setup: sqrt(P)*A * x = sqrt(P)*b, where sqrt(Pii) = sigma0/sigma(ii)
        Eigen::MatrixXd A = Eigen::MatrixXd(observations, parameters);
        Eigen::VectorXd b = Eigen::VectorXd(observations);
        Eigen::VectorXd x = Eigen::VectorXd(parameters);
        
        // Compute the mean epoch as mjd; all deltatimes are computed as differences
        // from this (mean) epoch.
        double mean_epoch { this->central_epoch().as_mjd() };

        // \todo would it be better to fill this column-wise??
        // Go ahead and form the A and b matrices.
        // We're gonna need some variables ...
        double dt, weight;
        std::size_t idx{0},     // index (i.e. row of A and b matrices)
                    counter{0}, // index of the current data point (m_data)
                    col{0};     // the column index
        epoch current_epoch;    // the current epoch

#ifdef DEBUG
        std::cout<<"\nNumber of data points: " << size();
        std::cout<<"\nNumber of epochs:      " << epochs();
        std::cout<<"\nA : "<<observations<<" * "<<parameters;
        std::cout<<"\nb : "<<observations;
        std::cout<<"\nStart: "<<first_epoch().stringify();
        std::cout<<" Stop: "<<last_epoch().stringify();
        std::cout<<" Mean: "<<central_epoch().stringify();
#endif

        for (auto it = m_data.cbegin(); it!= m_data.cend(); ++it)
        {
            // only include data that are not marked as 'skip'
            if ( !it->skip() ) {
                current_epoch = (*m_epochs)[counter];
                // delta days from central epoch as mjd.
                dt = current_epoch.as_mjd() - mean_epoch;
                // weight of observation
                weight = sigma0 / m_data[counter].sigma();
                //
                // Design Matrix A
                // ------------------------------------------------------------
                //
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
                // Set up jumps ...
                for (auto j = jumps.begin(); j != jumps.cend(); ++j) {
                    if ( *j >= current_epoch ) {
                        A(idx, col) = weight;
                    } else {
                        A(idx, col) = .0e0;
                    }
                    ++col;
                }
                // Set up velocity changes ...
                for (auto j = vel_changes.cbegin(); j != vel_changes.cend(); ++j) {
                    if ( *j >= current_epoch ) {
                        A(idx, col) = weight * (dt / 365.25);
                    } else {
                        A(idx, col) = .0e0;
                    }
                    ++col;
                }

                //
                // Observation Matrix (vector b)
                // ------------------------------------------------------------
                //
                b(idx) = m_data[counter].value() * weight;

#ifdef DEBUG
                assert(col == parameters);
#endif
                ++idx;
                col = 0;
            }
            ++counter;
        }
#ifdef DEBUG
        assert(counter <= size());
#endif
        
        // Solve via QR
        x = A.colPivHouseholderQr().solve(b);
        std::cout<<"\nLS solution vector: \n";
        for (std::size_t i=0;i<parameters;i++) printf("%+15.5f\n",x(i));

        // residual vector u = A*x - b
        Eigen::VectorXd u = Eigen::VectorXd(observations);
        u = A * x - b;
    
        // solution vector to std::vector
        std::vector<double> xvec;
        xvec.reserve(parameters);
        for (std::size_t ii = 0; ii < parameters; ii++) xvec.push_back(x(ii));

        // residuals as time-series
        timeseries<T, F> res {*this};
        idx = counter = 0;
        for (std::size_t i = 0; i < size(); i++) {
            if ( !m_data[i].skip() ) {
                data_point<F> dp { u(idx), m_data[i].sigma(),  m_data[i].flag() };
                ++idx;
                res[i] = dp;
            } else {
                res[i] = m_data[i];
            }
        }

        // should we print the model line to a file?
        /*
        if ( print_model_to ) {
            std::ofstream fout { (*print_model_to).c_str() };
            if ( !fout.is_open() ) {
                throw std::runtime_error(
                    "[ERROR] Failed to open file for writing: \""+(*print_model_to)+"\".");
            }
            std::vector<double> model_x, model_y;
            make_model_line(first_epoch(), last_epoch(), central_epoch(),
                xvec, jumps, vel_changes, periods, model_x, model_y);
        }
        */
        return xvec;
    }

    /// Return a const iterator to the first entry of the data points vector
    auto
    citer_start() const noexcept
    { return m_data.c_begin(); }
    
    /// Return a const iterator to the last+1 (=end) entry of the data points vector
    auto
    citer_stop() const noexcept
    { return m_data.c_end(); }

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

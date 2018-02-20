#ifndef __NGPT_TIMESERIES__
#define __NGPT_TIMESERIES__

// standard headers
#include <vector>
#include <iterator>
#include <algorithm>
#include <tuple>
#ifdef DEBUG
#include <iostream>
#include <iomanip>
#include <cstdio>
#endif

// Eigen headers
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/QR"
#include "eigen3/Eigen/Dense"

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "tsflag.hpp"
#include "model.hpp"

namespace ngpt
{

/// Forward declerations
template<class T, class F> class timeseries_iterator;
template<class T, class F> class timeseries_const_iterator;
template<class T, class F> class running_window;

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
template<class F>
    class data_point
{
public:

    /// Simplify the flag type.
    using tflag = ngpt::flag<F>;
    
    /// Constructor.
    explicit
    data_point(double val=0.0, double sigma=1e-3, tflag f=tflag{})
    noexcept
    : m_value{val},
      m_sigma{sigma},
      m_flag{f}
    {}

    /// (const) get the value.
    double
    value() const noexcept { return m_value; }

    /// get the value.
    double&
    value() noexcept { return m_value; }
    
    /// (const) get the sigma (std. dev).
    double
    sigma() const noexcept { return m_sigma; }
    
    /// get the sigma (std. dev).
    double&
    sigma() noexcept { return m_sigma; }
    
    /// (const) get the flag.
    tflag
    flag() const noexcept { return m_flag; }
    
    /// get the flag.
    tflag&
    flag() noexcept { return m_flag; }
    
    /// Should the data point be skipped/ignored ?
    /// @todo what the fuck is this? see also the detailed class description.
    bool
    skip() const noexcept { return ngpt::__skip__(this->m_flag); }

    /// equality operator
    bool
    operator==(const data_point& other) noexcept
    { return (m_value == other.m_value && m_sigma == other.m_sigma)
                                       && (m_flag == other.m_flag); }

private:
    double m_value; ///< The data point's value
    double m_sigma; ///< The data point's sigma (i.e. standard deviation)
    tflag  m_flag;  ///< The point's flag

}; // end class data_point


/// A generic time-series class.
/// Mean value and number of skipped points should always be correct
/// (i.e. updated).
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
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) of the time-series, will have a time-stamp of type
///           ngpt::datetime<T>.
/// @param F  An enumeration class to act as flag; each data point of the
///           time-series will have a flag of type ngpt::flag<F> (see class
///           data_point for details).
///
/// @example  test_ts.cpp
///
template<class T,
         class F,
         typename = std::enable_if_t<T::is_of_sec_type>
        >
    class timeseries
{
public:

    /// The specific datetime<T> class we will be using.
    using epoch = ngpt::datetime<T>;

    /// Simplify the flag type (usually consider ngpt::flag<pt_marker>). This
    /// acts like a flag for each time-series data point (a point can be an
    /// outlier or marked as skiped).
    using tflag = ngpt::flag<F>;

    /// The type of the time-series data points.
    using entry = ngpt::data_point<F>;

    /// The type of a time-series data entry, i.e. a pair of 
    /// <epoch, data_point>.
    using record = std::tuple<epoch&, entry&>;

    /// Assign the (non-skipped) time-series data_points and epochs to c-style
    /// double arrays. The input arrays, must be large enough to hold the data
    /// (they must have size >= data_pts() - skipped_pts()).
    /// The function also computes and assigns the time-series's data mean and
    /// variance.
    /// The average and variance are computed as:
    /// \f$\widetilde{x}_{n+1}=\frac{x_{n+1}+n*\widetilde{x_n}}{n+1} = \widetilde{x_n}+\frac{x_{n+1}-\widetilde{x_n}}{n+1}\f$
    /// \f$S_{n}=S_{n-1}+\frac{x_{n}-\widetilde{x}_{n-1}}{x_{n}-\widetilde{x_n}}\f$
    /// If we use weights, then the equations become:
    /// \f$\widetilde{x}_{n}=\widetilde{x}_{n-1}+\frac{w_n}{W_n}*(x_{n}-\widetilde{x}_{n-1})\f$
    /// \f$S_{n}=S_{n-1}+w_{n}(x_{n}-\widetilde{x}_{n-1})(x_{n}-\widetilde{x_n})\f$
    /// where \f$W_n=\sum_{i=1}^{n} w_i\f$
    std::size_t
    ts2array(double* epoch_arr, double* val_arr, double& ave, double& var)
    const
    {   
        std::size_t N {this->data_pts() - this->skipped_pts()};
        double prev_ave {0e0};
        auto   ts_start {this->cbegin()},
               ts_stop  {this->cend()};
        std::size_t index {0};
        double x,y;

        ave = var = 0e0;
        for (auto it = ts_start; it != ts_stop; ++it) {
            if ( !it.data().skip() ) {
                x = it.epoch().as_mjd();
                y = it.data().value();
                epoch_arr[index] = x;
                val_arr[index]   = y;
                prev_ave         = ave;
                ave += ((y-ave)/(index+1e0));
                var += ((y-prev_ave)*(y-ave));
                ++index;
            }
        }
        var = std::sqrt(var/N);

        return N;
    }
    
    /// Constructor. If a vector of epochs is passed in, then we know we have
    /// our epochs. In this case, reserve space (memory) for the data points.
    ///
    /// @param[in] epochs A pointer to a vector of ngpt::datetime<T> instances.
    ///
    /// @note Even though no data points are added to the time-series, enough
    ///       space is allocated to hold them (only allocated **not**
    ///       initialized!
    explicit
    timeseries(std::vector<epoch>* epochs=nullptr) noexcept
    : m_epochs(epochs),
      m_mean{0.0},
      m_skipped{0}
    {
        if ( m_epochs ) m_data.reserve(m_epochs->size());
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
    explicit
    timeseries(std::size_t size_hint) noexcept
    : m_epochs(nullptr),
      m_mean{0.0},
      m_skipped{0}
    {
        m_data.reserve(size_hint);
    }

    /// Get the (pointer to) epoch vector.
    std::vector<epoch>*&
    epoch_ptr() noexcept { return m_epochs; }

    /// Get the (pointer to) epoch vector (const version).
    const std::vector<epoch>*
    epoch_ptr() const noexcept { return m_epochs; }

    /// Get the data point (i.e. data_point<F>) at index i (const version).
    entry
    operator[](std::size_t i) const { return m_data[i]; }

    /// Get the data point (i.e. data_point<F>) at index i.
    entry&
    operator[](std::size_t i) { return m_data[i]; }
    
    /// Get the record (i.e. std::tuple<epoch&, entry&>) at index i.
    record&
    operator()(std::size_t i) { return std::tie((*m_epochs[i]), m_data[i]); }

    /// Get the epoch at a certain index (non-const)
    epoch&
    epoch_at(std::size_t i) { return m_epochs->operator[](i); }
    
    /// Get the epoch at a certain index (non-const)
    epoch
    epoch_at(std::size_t i) const { return m_epochs->operator[](i); }

    /// Get the mean value (average).
    /// @warning The mean value is **not** computed here; only the instance's
    ///          member m_value is returned.
    double
    mean() const noexcept { return m_mean; }

    /// Get the number of data points (all data points, regardless of their
    /// flag).
    std::size_t
    data_pts() const noexcept { return m_data.size(); }

    /// Get the number of epochs (all epochs, regardless if the corresponding
    /// data_point flags). If an epoch vector is not set, 0 is returned.
    std::size_t
    epochs() const noexcept { return m_epochs ? m_epochs->size() : 0; }
    
    /// Get the number of data points that are skipped. Not computation (check)
    /// is performed here; only the instance's member m_skipped is returned.
    /// This is the const version.
    std::size_t
    skipped_pts() const noexcept { return m_skipped; }
    
    /// Get the number of data points that are skipped. Not computation (check)
    /// is performed here; only the instance's member m_skipped is returned.
    std::size_t&
    skipped_pts() noexcept { return m_skipped; }

    /// Get the first epoch of the time-series (regardles of it's data_point
    /// flag).
    epoch
    first_epoch() const noexcept { return (*m_epochs)[0]; }
    
    /// Get the last epoch (regardles of it's data_point flag)
    epoch
    last_epoch() const noexcept { return (*m_epochs)[m_epochs->size()-1]; }

    /// Get the first epoch **NOT** skipped. Returns the value of the first,
    /// valid epoch and sets the idx parameter to its index.
    ///
    /// @param[out] idx  The index of the first not-skipped epoch.
    /// @todo            Should i allow this to throw?
    epoch
    first_valid_epoch(std::size_t& idx) const noexcept
    {
        auto it = std::find_if(std::cbegin(m_data), std::cend(m_data),
            [](const entry& i){return !(i.skip());});
        idx = std::distance(std::cbegin(m_data), it);
        return (*m_epochs)[idx];
    }

    /// Get the last epoch **NOT** skipped. Returns the value of the last,
    /// valid epoch and sets the idx parameter to its index.
    ///
    /// @param[out] idx  The index of the last not-skipped epoch.
    /// @todo            Should i allow this to throw?
    epoch
    last_valid_epoch(std::size_t& idx) const noexcept
    {
        auto it = std::find_if(std::crbegin(m_data), std::crend(m_data),
            [](const entry& i){return !(i.skip());});
        idx = std::distance(it, std::crend(m_data)) - 1;
        return (*m_epochs)[idx];
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
    timeseries
    (const timeseries& ts, std::size_t start=0, std::size_t end=0)
    : m_epochs(ts.m_epochs),
      m_mean{ts.m_mean},
      m_skipped{ts.m_skipped}
    {
        if (start || end) {
            if (start && !end) {
                end = ts.m_data.size();
            }
            if (end < start || end > ts.m_data.size()) {
                throw std::domain_error("timeseries: Invalid start/stop indexes for copy c'tor");
            }
            m_data.reserve(end-start);
            double sz;
            m_mean    = 0e0;
            m_skipped = 0;
            for (std::size_t i = start; i < end; ++i) {
                assert( i < ts.m_data.size() );
                m_data.emplace_back(ts.m_data[i]);
                if ( ts.m_data[i].skip() ) ++m_skipped;
                sz = i - start;
                m_mean = (ts[i].value() + sz*m_mean)/(sz+1e0);
            }
        } else {
            m_data = ts.m_data;
        }
    }

    /// Move constructor.
    /// @note The resluting time-serie's epoch vector (pointer), will be set to
    /// (point to) the original (i.e. ts). 
    timeseries(timeseries&& ts) noexcept
    : m_epochs(ts.m_epochs),
      m_mean{std::move(ts.m_mean)},
      m_data{std::move(ts.m_data)},
      m_skipped{std::move(ts.m_skipped)}
    {}

    /// Assignment operator.
    /// @warning Note that the epoch vector is not (deep) copied; it is only
    ///          set to point to the same epoch vector as the copied-from time-
    ///          series.
    timeseries&
    operator=(const timeseries& ts) noexcept
    {
        if (this != &ts) {
            m_epochs  = ts.m_epochs;
            m_mean    = ts.m_mean;
            m_data    = ts.m_data;
            m_skipped = ts.m_skipped;
        }
        return *this;
    }
    
    /// Move assignment operator.
    timeseries&
    operator=(timeseries&& ts) noexcept
    {
        if (this != &ts) {
            m_epochs = ts.m_epochs;
            m_mean   = std::move(ts.m_mean);
            m_data   = std::move(ts.m_data);
            m_skipped = std::move(ts.m_skipped);
        }
        return *this;
    }

    /// Split a time-series; return two new time-series in the interval:
    /// [0-idx) and [idx-end).
    ///
    /// @return A tuple (pair) containing the two new time-series.
    ///
    /// @todo   What the fuck should i do with the epochs of each sub-timeseries??
    auto
    split(std::size_t idx) const
    {
        timeseries left  (*this, 0, idx);
        timeseries right (*this, idx);
        return std::make_tuple(std::move(left), std::move(right));
    }

    /// Add a data point; returns the new mean value.
    /// @note   The instance's mean value is updated; so is the number of
    ///         skipped data points (if needed).
    /// @return The updated time-series mean value.
    double
    add_point(const entry& e)

    {
        double sz = static_cast<double>(m_data.size());
        m_data.emplace_back(e);
        if ( e.skip() ) {
            ++m_skipped;
        } else {
            m_mean = (e.value() + sz * m_mean)/(sz+1e0);
        }
        return m_mean;
    }

    /// Add a data point; returns the new mean value.
    /// @note   The instance's mean value is updated; so is the number of
    ///         skipped data points (if needed).
    /// @return The updated time-series mean value.
    double
    add_point(double val, double sigma=1e-3, tflag f=tflag{})
    {
        return this->add_point( entry{val, sigma, f} );
    }

    /// Mark a data point given its index.
    /// @param[in] index The index of the data point to be marked.
    /// @param[in] f     An enum (class instance) of type F. The enum f is
    ///                  added to the flag of the data point with index `index` 
    /// @note            The number of skipped data points is updated if needed.'
    /// @todo            WTF is up with __skip__ ??? And, somethin should be
    ///                  done with the epoch vectors
    void
    mark(std::size_t index, F f)/* noexcept noop, m_data[index] could throw */
    {
        tflag previous = m_data[index].flag();
        m_data[index].flag().set(f);
        if ( !__skip__(previous) && __skip__(tflag{f}) ) ++m_skipped;
    }

    /// Compute the mean (i.e. central epoch). This version uses the very 
    /// first and last epochs to compute the mean, regardless if they are 
    /// marked as unused. The mean epoch is obviously half the distance
    /// between the first and last epochs.
    epoch
    central_epoch() const noexcept
    {
        auto delta_dt  = ngpt::delta_date(last_epoch(), first_epoch());
        auto central_epoch { first_epoch() };
        central_epoch += (delta_dt / 2);
        return central_epoch;
    }
    
    /// Compute the mean (i.e. central epoch). This version uses the first and
    /// last epochs that are not marked as unused, to compute the mean epoch. 
    /// I.e., if the first 10 epochs --data points-- are marked as outliers, 
    /// then they shall not be used to compute the mean.
    /// The mean epoch is obviously half the distance between the first and 
    /// last epochs.
    epoch
    central_valid_epoch() const noexcept
    {
        auto start_epoch = this->first_valid_epoch();
        auto delta_dt    = ngpt::delta_date(start_epoch, last_valid_epoch());
        start_epoch     += delta_dt;
        return start_epoch;
    }

    /// @todo this should be const, but there is a problem in line 564
    /*
    auto
    detrend(double& x0, double& vx, bool mark_outliers=false, double sigma0=1e-3) 
    {
        // Construct a linear model
        ngpt::ts_model<T> model;
        std::size_t parameters = model.parameters();
        std::size_t observations = m_data.size() - m_skipped;

        // Set up the matrices; the model is: A * x = b
        // since we have uncorrelated variance matric though, we will actually
        // setup: sqrt(P)*A * x = sqrt(P)*b, where sqrt(Pii) = sigma0/sigma(ii)
        Eigen::MatrixXd A = Eigen::MatrixXd(observations, parameters);
        Eigen::VectorXd b = Eigen::VectorXd(observations);
        Eigen::VectorXd x = Eigen::VectorXd(parameters);
        
        // Compute the mean epoch as mjd; all deltatimes are computed as differences
        // from this (mean) epoch.
        model.mean_epoch() = this->central_epoch();
    
        double dt, weight;
        std::size_t idx{0},       // index (i.e. row of A and b matrices)
                    counter{0};   // index of the current data point (m_data)
        datetime<T> current_epoch;// the current epoch
        timeseries_const_iterator<T, F> it = cbegin(),
                                    it_end = cend();
        data_point<F> entry;

        for (; it != it_end; ++it) {
            entry         = it.data();
            current_epoch = it.epoch();
            // only include data that are not marked as 'skip'
            if ( !entry.skip() ) {
                dt = current_epoch.as_mjd() - model.mean_epoch().as_mjd();
                weight = sigma0 / entry.sigma();
                model.assign_row(A, b, current_epoch, weight, dt, entry.value(), idx);
                ++idx;
            }
            ++counter;
        }
        assert( counter == data_pts() );

        // Solve via QR
        x = A.colPivHouseholderQr().solve(b);

        // Variance-Covariance Matrix
        int N = (A.transpose()*A).rows();
        auto Q = (A.transpose() * A).fullPivLu().solve(Eigen::MatrixXd::Identity(N,N));
        std::cout<<"\n\tVar-Covar Info:";
        std::cout<<"\n\tQ is of type:"<<Q.rows()<<"x"<<Q.cols();

        // residual vector u = A*x - b; note that the residual vector may not
        // have the same size as the (original) time-series. Instead, it has
        // a size of: original_ts.size() - original_ts.skipped_pts().
        //
        // IMPORTANT
        // ---------------------------
        // Note that at this point the residuals are scaled according to each
        // data points weight! (We actualy solved not Ax=b but QAx=Qb)
        Eigen::VectorXd u = Eigen::VectorXd(observations);
        u = A * x - b;
    
        // assign solution vector to the model.
        model.assign_solution_vector(x);

        // Cast residuals to time-series and compute a-posteriori std. dev
        // The resulting residual ts will have the same size as the original ts,
        // where data points marked as 'skipped' will have their original value.
        double post_std_dev {0e0};
        timeseries<T, F> res {*this};
        res.epoch_ptr() = this->epoch_ptr();
        idx = counter = 0;
        double residual;
        for (std::size_t i = 0; i < epochs(); i++) {
            if ( !m_data[i].skip() ) {
                residual = u(idx)/(sigma0/m_data[i].sigma());
                data_point<F> dp { residual, m_data[i].sigma(),  m_data[i].flag() };
                post_std_dev += residual*residual;
                res[i] = dp;
                ++idx;
            } else {
                res[i] = m_data[i];
            }
        }
        post_std_dev = std::sqrt(post_std_dev)
            /(double)(idx-parameters);

        x0 = model.x0();
        vx = model.vx();
        std::cout<<"\n[DEBUG] De-trended component, with velocity: "<<vx;
        if (mark_outliers) {
            datetime_interval<T> window {modified_julian_day{90}, T{0}};
            nikolaidis(res, *this, window);
        }
        return res;
    }*/

    /// Given a model, ...., he... well solve for it!
    /// @todo document a little, just a bit, better.
    /// @param[in] model           An instance of type md_model; the parameter
    ///                            estimation process, is based on this model.
    /// @param[out] post_std_dev   A-posteriori std. deviation.
    /// @param[in] sigma0          A-priori std. deviation.
    /// @param[in] set_model_epoch A boolean variable, to denote:
    ///                            set_model_epoch | Action
    ///                            ----------------|--------------------------
    ///                             true           | To estimate the paramters,
    ///                                            | the timeserie's mean epoch
    ///                                            | is used; this value is also
    ///                                            | assigned to the model's
    ///                                            | mean_epoch member variable.
    ///                             ---------------|---------------------------
    ///                              false         | To estimate the paramters,
    ///                                            | the model's mean_epoch
    ///                                            | is used (as central epoch
    ///                                            | for the computations).
    ///
    /// @note
    ///        - Number of model parametrs must be positive (>0)
    ///        - Number of observations must be larger than number of
    ///          parameters
    ///        - Only non-skiiped data points are considered (i.e. all
    ///          data-points for which the skip() function returns false value.
    auto
    qr_ls_solve(
        ngpt::ts_model<T>&    model,
        double&               post_std_dev,
        double sigma0         = 1e-03,
        bool   set_model_epoch= false, /* set or not the mean epoch for the (input) model */
        bool   mark_outliers  = true
    )
    {
        // Number of parameters to be estimated:
        std::size_t parameters = model.parameters();
        assert( parameters > 0 );

        // Number of observations (ommiting the ones to be skipped)
        std::size_t observations = m_data.size() - m_skipped;
        assert( observations > parameters );
        
        // Set up the matrices; the model is: A * x = b
        // since we have uncorrelated variance matric though, we will actually
        // setup: sqrt(P)*A * x = sqrt(P)*b, where sqrt(Pii) = sigma0/sigma(ii)
        Eigen::MatrixXd A = Eigen::MatrixXd(observations, parameters);
        Eigen::VectorXd b = Eigen::VectorXd(observations);
        Eigen::VectorXd x = Eigen::VectorXd(parameters);
        
        // Compute the mean epoch as mjd; all deltatimes are computed as differences
        // from this (mean) epoch.
        double mean_epoch { this->central_epoch().as_mjd() };
        if ( set_model_epoch ) {
            model.mean_epoch() = this->central_epoch();
        } else {
            mean_epoch = model.mean_epoch().as_mjd();
        }
    
        double dt, weight;
        std::size_t idx{0},       // index (i.e. row of A and b matrices)
                    counter{0};   // index of the current data point (m_data)
        datetime<T> current_epoch;// the current epoch
        timeseries_const_iterator<T, F> it = cbegin(),
                                    it_end = cend();
        data_point<F> entry;

        for (; it != it_end; ++it) {
            entry         = it.data();
            current_epoch = it.epoch();
            // only include data that are not marked as 'skip'
            if ( !entry.skip() ) {
                dt = current_epoch.as_mjd() - mean_epoch;
                weight = sigma0 / entry.sigma();
                model.assign_row(A, b, current_epoch, weight, dt, entry.value(), idx);
                ++idx;
            }
            ++counter;
        }
        assert( counter == data_pts() );

        // Solve via QR
        // x = A.colPivHouseholderQr().solve(b);
        // Solve via SVD
        x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

        // residual vector u = A*x - b; note that the residual vector may not
        // have the same size as the (original) time-series. Instead, it has
        // a size of: original_ts.size() - original_ts.skipped_pts().
        //
        // IMPORTANT
        // ---------------------------
        // Note that at this point the residuals are scaled according to each
        // data points weight! (We actualy solved not Ax=b but QAx=Qb)
        Eigen::VectorXd u = Eigen::VectorXd(observations);
        u = A * x - b;
    
        // assign solution vector to the model.
        model.assign_solution_vector(x);

        // Cast residuals to time-series and compute a-posteriori std. dev
        // The resulting residual ts will have the same size as the original ts,
        // where data points marked as 'skipped' will have their original value.
        post_std_dev = 0;
        timeseries<T, F> res {*this};
        res.epoch_ptr() = this->epoch_ptr();
        idx = counter = 0;
        double residual;
        for (std::size_t i = 0; i < epochs(); i++) {
            if ( !m_data[i].skip() ) {
                residual = u(idx)/(sigma0/m_data[i].sigma()); /* is this correct, or should it be u(idx)/(sigma0^2/m_data[i].sigma()^2) */
                data_point<F> dp { residual, m_data[i].sigma(),  m_data[i].flag() };
                post_std_dev += residual*residual;
                res[i] = dp;
                ++idx;
            } else {
                res[i] = m_data[i];
            }
        }
        post_std_dev = std::sqrt(post_std_dev)
            /(double)(idx-parameters);
        
        // --------------------------------------------------------------------
        // Variance-Covariance Matrix
        // see https://forum.kde.org/viewtopic.php?f=74&t=139476
        // --------------------------------------------------------------------
        int N = (A.transpose()*A).rows();
        Eigen::MatrixXd Q = Eigen::MatrixXd(parameters, parameters);
        // Q = (A.transpose() * A).fullPivLu().solve(Eigen::MatrixXd::Identity(N,N));
        Q = (A.transpose() * A).ldlt().solve(Eigen::MatrixXd::Identity(N,N));
        
        // std::cout<<"\n\tVar-Covar Info:";
        // std::cout<<"\n\tQ is of type:"<<Q.rows()<<"x"<<Q.cols()<<" (N="<<N<<")";
        // for (int i=0; i<N; i++) {
        //    std::cout<<"\n\t"<<x(i)<<" +/- "<<post_std_dev*std::sqrt( Q(i,i) ) << " ratio: "<<x(i) / (post_std_dev*std::sqrt( Q(i,i) ));
        // }

        // apply outlier detection algorithm and mark them
        if (mark_outliers) {
            datetime_interval<T> window {modified_julian_day{90}, T{0}};
            nikolaidis(res, *this, window);
        }
        return res;
    }

    /// Return a const iterator to the first entry of the data points vector
    /// (i.e. return m_data.cbegin())
    typename std::vector<entry>::const_iterator
    citer_start() const noexcept
    { return m_data.cbegin(); }

    /// Return an (non-const) iterator to the first entry of the data points
    /// vector (i.e. return m_data.begin())
    typename std::vector<entry>::iterator
    iter_start() noexcept
    { return m_data.begin(); }
    
    /// Return a const iterator to the last+1 (=end) entry of the data points
    /// vector (i.e. return m_data.cend())
    typename std::vector<entry>::const_iterator
    citer_stop() const noexcept
    { return m_data.cend(); }
    
    /// Return a (non-const) iterator to the last+1 (=end) entry of the data
    /// points vector (i.e. return m_data.end())
    typename std::vector<entry>::iterator
    iter_stop() noexcept
    { return m_data.end(); }

    /// Return a timeseries_iterator to the first element of this instance
    /// (reagrdless of its flag).
    /// @todo is it noexcept?
    timeseries_iterator<T, F>
    begin()
    { return timeseries_iterator<T, F> {*this}; }

    /// Return a timeseries_iterator to the one-past-the-end element of this 
    /// instance (reagrdless of its flag).
    /// @todo is it noexcept?
    timeseries_iterator<T, F>
    end()
    {
        timeseries_iterator<T, F> it {*this};
        it.end();
        return it;
    }
    
    /// Return a cosnt timeseries_iterator to the first element of this instance
    /// (reagrdless of its flag).
    /// @todo is it noexcept?
    timeseries_const_iterator<T, F>
    cbegin() const
    { return timeseries_const_iterator<T, F> {*this}; }

    /// Return a const timeseries_iterator to the one-past-the-end element of 
    /// this instance (reagrdless of its flag).
    /// @todo is it noexcept?
    timeseries_const_iterator<T, F>
    cend() const
    {
        timeseries_const_iterator<T, F> it {*this};
        it.end();
        return it;
    }

    /// Return an iterator pointing to the first element in the time-series 
    /// with a time-stamp that is greater than the given value, or last if no 
    /// such element is found. 
    /// @param[in] t    The input time-stamp; we are searching for the first
    ///                 time-series element with time-stamp > t.
    /// @param[out] idx The index of the time-series element with time-stamp > t,
    ///                 If no such element is found, idx is set equal to the
    ///                 total number of epochs.
    /// @return         A timeseries_const_iterator pointing to the the first
    ///                 element in the time-series with a time-stamp > t; if not
    ///                 found, a const iterator to end.
    timeseries_const_iterator<T, F>
    upper_bound(const datetime<T>& t, std::size_t& idx)
    const noexcept
    {
        auto it = std::upper_bound(std::cbegin(*m_epochs),
                                   std::cend(*m_epochs),
                                   t);
        if ( it == m_epochs->cend() ) {
            idx = m_epochs->size();
            return this->cend();
        } else {
            idx = std::distance(std::cbegin(*m_epochs), it);
            timeseries_const_iterator<T,F> itt (*this);
            return ( itt+idx );
        }
    }
    
    /// Return an iterator pointing to the first element in the time-series 
    /// with a time-stamp that is not less than the given value, or last if no 
    /// such element is found. 
    /// @param[in] t    The input time-stamp; we are searching for the first
    ///                 time-series element with time-stamp >= t.
    /// @param[out] idx The index of the time-series element with time-stamp >= t,
    ///                 If no such element is found, idx is set equal to the
    ///                 total number of epochs.
    /// @return         A timeseries_const_iterator pointing to the the first
    ///                 element in the time-series with a time-stamp >= t; if not
    ///                 found, a const iterator to end.
    timeseries_const_iterator<T, F>
    lower_bound(const datetime<T>& t, std::size_t& idx)
    const noexcept
    {
        auto it = std::lower_bound(std::cbegin(*m_epochs),
                                   std::cend(*m_epochs),
                                   t);
        if ( it == m_epochs->cend() ) {
            idx = m_epochs->size();
            return this->cend();
        } else {
            idx = std::distance(std::cbegin(*m_epochs), it);
            timeseries_const_iterator<T,F> itt (*this);
            return ( itt+idx );
        }
    }

    running_window<T, F>
    rw_begin(const datetime_interval<T>& window)
    {
        running_window<T, F> it {*this, window};
        it.begin();
        return it;
    }

    running_window<T, F>
    rw_end(const datetime_interval<T>& window)
    {
        running_window<T, F> it {*this, window};
        it.end();
        return it;
    }

#ifdef DEBUG
    /// Write the time-series to an output stream. This is only very basic, it
    /// needs work.
    /// @todo 
    ///       - make this function an '<<' operator.
    ///       - this means that the flag has an << operator
    std::ostream&
    dump(std::ostream& os) const
    {
        auto iter = this->cbegin();

        for (; iter != this->cend(); ++iter) {
            os <<  iter.epoch().as_mjd() << " "
                << iter.data().value() << " " 
                << iter.data().sigma() << " " 
                << iter.data().flag()  << "\n";
        }
        return os;
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
    std::size_t m_skipped;

}; // class timeseries

/// The time_series class, should be iterable, i.e. the user should be able to
/// 'walk through' an instance, data-point by data-point. This class, enables
/// this feature.
/// It is an effort to make something like an std::vector<T>::(const)_)iterator.
/// The data-points flags do not play any role in the iteration, i.e. when we
/// iterate we walk through all data-points regardless of their flags.
/// The template parameters (T and F), should be the same as for the actual
/// timeseries we want to iterate.
/// The internals of this class are pretty simple; it just holds a (reference
/// to) time_series instance and two iterators, one for the epoch vector and
/// one for the data_point vector (both are of type std::vector).
///
/// @todo Do i really need the epoch iterator as a member ?? Maybe just use
///       pointer arithmetic for this! This indeed can be done (see the code
///       below that is commented out), but shows no significant difference
///       in efficiency (it is not faster!).
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) of the time-series, will have a time-stamp of type
///           ngpt::datetime<T>.
/// @param F  An enumeration class to act as flag; each data point of the
///           time-series will have a flag of type ngpt::flag<F> (see class
///           data_point for details).
///
/// @example  test_ts.cpp
///
template<class T, class F>
    class timeseries_iterator
{
public:

    /// The specific datetime<T> class we will be using.
    using epoch_td    = ngpt::datetime<T>;
    using interval_td = ngpt::datetime_interval<T>;
    using tflag       = ngpt::flag<F>;
    using entry       = ngpt::data_point<F>;
    using record      = std::tuple<epoch_td&, entry&>;

    /// Constructor given a time-series instance.
    explicit
    timeseries_iterator(timeseries<T, F, typename std::enable_if_t<T::is_of_sec_type>>& ts)
    : m_timeseries{ts},
      m_data_iter{ts.iter_start()},
      m_epoch_iter{ts.epoch_ptr()->begin()}
    { assert( ts.epoch_ptr() && ts.data_pts() == ts.epochs() ); }

    /// Copy constructor
    timeseries_iterator(const timeseries_iterator& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {}

    /// Move constructor
    timeseries_iterator(timeseries_iterator&& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {}

    /// Assignment operator
    timeseries_iterator&
    operator=(const timeseries_iterator& it) noexcept
    {
        if (*this != it) {
            m_timeseries = it.m_timeseries;
            m_data_iter  = it.m_data_iter;
            m_epoch_iter = it.m_epoch_iter;
        }
        return *this;
    }
    
    /// Move assignment operator
    timeseries_iterator&
    operator=(timeseries_iterator&& it) noexcept
    {
        if (*this != it) {
            m_timeseries = std::move(it.m_timeseries);
            m_data_iter  = std::move(it.m_data_iter);
            m_epoch_iter = std::move(it.m_epoch_iter);
        }
        return *this;
    }

    /// Set the iterator to the beggining of the time-series
    void
    begin() noexcept
    {
        m_data_iter  = m_timeseries.iter_start();
        m_epoch_iter = m_timeseries.epoch_ptr()->begin();
    }

    /// Set the iterator to one-past-the-end of the time-series
    void
    end() noexcept
    {
        m_data_iter  = m_timeseries.iter_stop();
        m_epoch_iter = m_timeseries.epoch_ptr()->end();
    }

    /// Move the iterator to the next data_point/epoch.
    void
    advance() noexcept
    {
        ++m_data_iter;
        ++m_epoch_iter;
    }

    /// Return the index of the iterator (i.e. the number of the current
    /// epoch/data-point)
    std::size_t
    index() noexcept
    { return std::distance(m_timeseries.iter_start(), m_data_iter); }

    /// Equality operator (i.e. both iterators point to the same element)
    bool
    operator==(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter  == it.m_data_iter
              && m_epoch_iter == it.m_epoch_iter );
    }

    /// InEquality operator (i.e. iterators do not point to the same element)
    bool
    operator!=(const timeseries_iterator& it) const noexcept
    { return !((*this) == it); }

    /// Bigger-than operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator>(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter  > it.m_data_iter
              && m_epoch_iter > it.m_epoch_iter );
    }
    
    /// Bigger-or-equal-to operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator>=(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter  >= it.m_data_iter
              && m_epoch_iter >= it.m_epoch_iter );
    }
    
    /// Less than operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator<(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter  < it.m_data_iter
              && m_epoch_iter < it.m_epoch_iter );
    }
    
    /// Less-or-equal-to operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator<=(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter  <= it.m_data_iter
              && m_epoch_iter <= it.m_epoch_iter );
    }

    /// Move the iterator n places (indexes) forward.
    timeseries_iterator // pointer arithmetic
    operator+(int n) const noexcept
    {
        timeseries_iterator ti {*this};
        ti.m_data_iter += n;
        ti.m_epoch_iter += n;
        return ti;
    }
    /// Move the iterator n places (indexes) backwards.
    timeseries_iterator // pointer arithmetic
    operator-(int n) const noexcept
    {
        timeseries_iterator ti {*this};
        ti.m_data_iter -= n;
        ti.m_epoch_iter -= n;
        return ti;
    }

    // Prefix operator (addition).
    timeseries_iterator&
    operator++()
    {
        ++m_data_iter;
        ++m_epoch_iter;
        return *this;
    }
    
    // Prefix operator (subtraction)
    timeseries_iterator& operator--()
    {
        --m_data_iter;
        --m_epoch_iter;
        return *this;
    }

    // Postfix operator (addition)
    timeseries_iterator
    operator++(int)
    {
        auto tmp {*this};
        this->operator++();
        return tmp;
    }

    /// Return the distance (i.e. index difference) between two iterators.
    /// Obviously, they must belong to the same time-series.
    int
    distance_from(const timeseries_iterator& it) const
    {
        assert( &m_timeseries == &it.m_timeseries );
        return std::distance(it.m_data_iter, this->m_data_iter);
    }

    /// Return the time-delta between two iterators.
    interval_td
    delta_time(const timeseries_iterator& it) const
    {
        assert( &m_timeseries == &(it.m_timeseries) );
        return m_epoch_iter->delta_date( *(it.m_epoch_iter) );
    }

    /// Return the epoch the iterator instance is currently at.
    epoch_td&
    epoch() noexcept { return *m_epoch_iter;}

    /// Return the entry the iterator instance is currently at.
    entry&
    data() noexcept { return *m_data_iter; }
    
    /// Mark the data point this instance is currently at.
    void
    mark(F f)
    {
        tflag previous = m_data_iter->flag();
        m_data_iter->flag().set(f);
        if ( !__skip__(previous) && __skip__(tflag{f}) ) {
            m_timeseries.skipped_pts() = m_timeseries.skipped_pts() + 1;
        }
    }

    /// Return (a reference to) the time_series instance this iterator points to.
    timeseries<T, F>&
    timeseries_ref() noexcept
    { return m_timeseries; }

private:
    /// The (reference to) time_series instance this iterator points to.
    timeseries<T, F>&                        m_timeseries;
    /// Iterator (i.e std::vector<>::iterator) of the data-points vector
    typename std::vector<entry>::iterator    m_data_iter;
    /// Iterator (i.e std::vector<>::iterator) of the epochs vector
    typename std::vector<epoch_td>::iterator m_epoch_iter;
}; // class timeseries_iterator

/*
/// The time_series class, should be iterable, i.e. the user should be able to
/// 'walk through' an instance, data-point by data-point. This class, enables
/// this feature.
/// It is an effort to make something like an std::vector<T>::(const)_)iterator.
/// The data-points flags do not play any role in the iteration, i.e. when we
/// iterate we walk through all data-points regardless of their flags.
/// The template parameters (T and F), should be the same as for the actual
/// timeseries we want to iterate.
/// The internals of this class are pretty simple; it just holds a (reference
/// to) time_series instance and one index (for both the poch and the data
/// points vector).
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) of the time-series, will have a time-stamp of type
///           ngpt::datetime<T>.
/// @param F  An enumeration class to act as flag; each data point of the
///           time-series will have a flag of type ngpt::flag<F> (see class
///           data_point for details).
///
/// @example  test_ts.cpp
template<class T, class F>
    class timeseries_iterator
{
public:

    /// The specific datetime<T> class we will be using.
    using epoch_td    = ngpt::datetime<T>;
    using interval_td = ngpt::datetime_interval<T>;
    using tflag       = ngpt::flag<F>;
    using entry       = ngpt::data_point<F>;
    using record      = std::tuple<epoch_td&, entry&>;

    /// Constructor given a time-series instance.
    explicit
    timeseries_iterator(timeseries<T, F, typename std::enable_if_t<T::is_of_sec_type>>& ts)
    : m_timeseries{ts},
      m_idx{0}
    { assert( ts.epoch_ptr() && ts.data_pts() == ts.epochs() ); }

    /// Copy constructor
    timeseries_iterator(const timeseries_iterator& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_idx{it.m_idx}
    {}

    /// Move constructor
    timeseries_iterator(timeseries_iterator&& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_idx{it.m_idx}
    {}

    /// Assignment operator
    timeseries_iterator&
    operator=(const timeseries_iterator& it) noexcept
    {
        if (*this != it) {
            m_timeseries = it.m_timeseries;
            m_idx        = it.m_idx;
        }
        return *this;
    }
    
    /// Move assignment operator
    timeseries_iterator&
    operator=(timeseries_iterator&& it) noexcept
    {
        if (*this != it) {
            m_timeseries = std::move(it.m_timeseries);
            m_idx        = std::move(it.m_idx);
        }
        return *this;
    }

    /// Set the iterator to the beggining of the time-series
    void
    begin() noexcept
    { m_idx = 0; }

    /// Set the iterator to one-past-the-end of the time-series
    void
    end() noexcept
    { m_idx = m_timeseries.data_pts(); }

    /// Move the iterator to the next data_point/epoch.
    void
    advance() noexcept { ++m_idx; }

    /// Return the index of the iterator (i.e. the number of the current
    /// epoch/data-point)
    std::size_t
    index() noexcept { return m_idx; }

    /// Equality operator (i.e. both iterators point to the same element)
    bool
    operator==(const timeseries_iterator& it) const noexcept
    { return m_idx == it.m_idx; }

    /// InEquality operator (i.e. iterators do not point to the same element)
    bool
    operator!=(const timeseries_iterator& it) const noexcept
    { return !((*this) == it); }

    /// Bigger-than operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator>(const timeseries_iterator& it) const noexcept
    { return m_idx > it.m_idx; }
    
    /// Bigger-or-equal-to operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator>=(const timeseries_iterator& it) const noexcept
    { return m_idx >= it.m_idx; }
    
    /// Less than operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator<(const timeseries_iterator& it) const noexcept
    { return m_idx < it.m_idx; }
    
    /// Less-or-equal-to operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator<=(const timeseries_iterator& it) const noexcept
    { return m_idx <= it.m_idx; }

    /// Move the iterator n places (indexes) forward.
    timeseries_iterator // pointer arithmetic
    operator+(int n) const noexcept
    {
        timeseries_iterator ti {*this};
        ti.m_idx += n;
        return ti;
    }
    /// Move the iterator n places (indexes) backwards.
    timeseries_iterator // pointer arithmetic
    operator-(int n) const noexcept
    {
        timeseries_iterator ti {*this};
        ti.m_idx -= n;
        return ti;
    }

    // Prefix operator (addition).
    timeseries_iterator&
    operator++()
    { 
        ++m_idx;
        return *this;
    }
    
    // Prefix operator (subtraction)
    timeseries_iterator& operator--()
    {
        --m_idx;
        return *this;
    }

    // Postfix operator (addition)
    timeseries_iterator
    operator++(int)
    {
        auto tmp {*this};
        this->operator++();
        return tmp;
    }

    /// Return the distance (i.e. index difference) between two iterators.
    /// Obviously, they must belong to the same time-series.
    int
    distance_from(const timeseries_iterator& it) const
    {
        assert( &m_timeseries == &it.m_timeseries );
        return m_idx - it.m_idx;
    }

    /// Return the time-delta between two iterators.
    interval_td
    delta_time(const timeseries_iterator& it) const
    {
        assert( &m_timeseries == &(it.m_timeseries) );
        auto t1 { m_timeseries.epoch_at(m_idx) };
        auto t2 { m_timeseries.epoch_at(it.m_idx) };
        return t1.delta_date( t2 );
    }

    /// Return the epoch the iterator instance is currently at.
    epoch_td&
    epoch() noexcept { return  m_timeseries.epoch_ptr[m_idx]; }

    /// Return the entry the iterator instance is currently at.
    /// @warning Via tis function, you can have access to the data_point at
    ///          any given (valid) index. However, doing this, will **not**
    ///          automatically update the internals of the time-series. E.g.
    ///          you could use this function to mark a data_point as outlier,
    ///          but the timeseries's m_skipped (i.e. the number of skipped
    ///          data points) member variable will **not** be updated!
    ///
    /// @todo    See the warning above; i really need to fix this some way. Also,
    ///          document this behaviour somewhere where it can be read ALWAYS.
    entry&
    data() noexcept { return m_timeseries[m_idx]; }
    
    /// Mark the data point this instance is currently at.
    void
    mark(F f)
    {
        tflag previous = m_timeseries[m_idx].flag();
        m_timeseries[m_idx].flag().set(f);
        if ( !__skip__(previous) && __skip__(tflag{f}) ) {
            m_timeseries.skipped_pts() = m_timeseries.skipped_pts() + 1;
        }
    }

    /// Return (a reference to) the time_series instance this iterator points to.
    timeseries<T, F>&
    timeseries_ref() noexcept
    { return m_timeseries; }

private:
    /// The (reference to) time_series instance this iterator points to.
    timeseries<T, F>&                        m_timeseries;
    long                                     m_idx;
}; // class timeseries_iterator
*/

/// The time_series class, should be iterable, i.e. the user should be able to
/// 'walk through' an instance, data-point by data-point. This class, enables
/// this feature.
/// It is an effort to make something like an std::vector<T>::const_iterator.
/// The data-points flags do not play any role in the iteration, i.e. when we
/// iterate we walk through all data-points regardless of their flags.
/// The template parameters (T and F), should be the same as for the actual
/// timeseries we want to iterate.
/// Remember, this a **const** iterator; you can't use an instance of this
/// class to alter the values of the timeseries (either the time or data_points).
/// The internals of this class are pretty simple; it just holds a (reference
/// to) time_series instance and two iterators, one for the epoch vector and
/// one for the data_point vector (both are of type std::vector).
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) of the time-series, will have a time-stamp of type
///           ngpt::datetime<T>.
/// @param F  An enumeration class to act as flag; each data point of the
///           time-series will have a flag of type ngpt::flag<F> (see class
///           data_point for details).
///
/// @example test_ts.cpp
///
template<class T, class F>
    class timeseries_const_iterator
{
public:

    /// The specific datetime<T> class we will be using.
    using epoch_td    = ngpt::datetime<T>;
    using interval_td = ngpt::datetime_interval<T>;
    using tflag       = ngpt::flag<F>;
    using entry       = ngpt::data_point<F>;
    using record      = std::tuple<epoch_td&, entry&>;

    /// Constructor
    explicit
    timeseries_const_iterator(const timeseries<T, F, 
        typename std::enable_if_t<T::is_of_sec_type>>& ts)
    : m_timeseries{ts},
      m_data_iter{ts.citer_start()},
      m_epoch_iter{ts.epoch_ptr()->cbegin()}
    { assert( ts.epoch_ptr() && (ts.data_pts() == ts.epochs()) ); }

    /// Copy constructor
    timeseries_const_iterator(const timeseries_const_iterator& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {}

    /// Move constructor
    timeseries_const_iterator(timeseries_const_iterator&& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {}

    /// Assignment operator
    timeseries_const_iterator&
    operator=(const timeseries_const_iterator& it) noexcept
    {
        if (*this != it) {
            m_timeseries = it.m_timeseries;
            m_data_iter  = it.m_data_iter;
            m_epoch_iter = it.m_epoch_iter;
        }
        return *this;
    }
    
    /// Move assignment operator
    timeseries_const_iterator&
    operator=(timeseries_const_iterator&& it) noexcept
    {
        if (*this != it) {
            m_timeseries = std::move(it.m_timeseries);
            m_data_iter  = std::move(it.m_data_iter);
            m_epoch_iter = std::move(it.m_epoch_iter);
        }
        return *this;
    }

    /// Set the iterator to the beggining of the time-series
    void
    begin() noexcept
    {
        m_data_iter  = m_timeseries.citer_start();
        m_epoch_iter = m_timeseries.epoch_ptr()->cbegin();
    }

    /// Set the iterator to one-past-the-end of the time-series
    void
    end() noexcept
    {
        m_data_iter  = m_timeseries.citer_stop();
        m_epoch_iter = m_timeseries.epoch_ptr()->cend();
    }

    /// Move the iterator to the next data_point/epoch.
    void
    advance() noexcept
    {
        ++m_data_iter;
        ++m_epoch_iter;
    }

    /// Return the index of the iterator (i.e. the number of the current
    /// epoch/data-point)
    std::size_t
    index() noexcept
    { return std::distance(m_timeseries.citer_start(), m_data_iter); }

    /// Equality operator (i.e. both iterators point to the same element)
    bool
    operator==(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter  == it.m_data_iter
              && m_epoch_iter == it.m_epoch_iter );
    }
    
    /// InEquality operator (i.e. iterators do not point to the same element)
    bool
    operator!=(const timeseries_const_iterator& it) const noexcept
    { return !((*this) == it); }

    /// Bigger-than operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator>(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter  > it.m_data_iter
              && m_epoch_iter > it.m_epoch_iter );
    }
    
    /// Bigger-or-equal-to operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator>=(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter  >= it.m_data_iter
              && m_epoch_iter >= it.m_epoch_iter );
    }
    
    /// Less than operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator<(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter < it.m_data_iter
              && m_epoch_iter < it.m_epoch_iter );
    }
    
    /// Less-or-equal-to operator. Checks the order (i.e. indexes) of the two
    /// iterators.
    bool
    operator<=(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter  <= it.m_data_iter
              && m_epoch_iter <= it.m_epoch_iter );
    }

    /// Move the iterator n places (indexes) forward.
    timeseries_const_iterator // pointer arithmetic
    operator+(int n) const noexcept
    {
        timeseries_const_iterator ti {*this};
        ti.m_data_iter  += n;
        ti.m_epoch_iter += n;
        return ti;
    }

    /// Move the iterator n places (indexes) backwards.
    timeseries_const_iterator // pointer arithmetic
    operator-(int n) const noexcept
    {
        timeseries_const_iterator ti {*this};
        ti.m_data_iter  -= n;
        ti.m_epoch_iter -= n;
        return ti;
    }

    // Prefix operator
    timeseries_const_iterator&
    operator++()
    {
        ++m_data_iter;
        ++m_epoch_iter;
        return *this;
    }
    
    // Prefix operator
    timeseries_const_iterator&
    operator--()
    {
        --m_data_iter;
        --m_epoch_iter;
        return *this;
    }

    // Postfix operator
    timeseries_const_iterator
    operator++(int)
    {
        auto tmp {*this};
        this->operator++();
        return tmp;
    }

    /// Return the distance (i.e. index difference) between two iterators.
    /// Obviously, they must belong to the same time-series.
    int
    distance_from(const timeseries_const_iterator& it) const
    {
        assert( &m_timeseries == &it.m_timeseries );
        return std::distance(it.m_data_iter, this->m_data_iter);
    }

    /// Return the time-delta between two iterators.
    interval_td
    delta_time(const timeseries_const_iterator& it) const
    {
        assert( &m_timeseries == &(it.m_timeseries) );
        return m_epoch_iter->delta_date( *(it.m_epoch_iter) );
    }

    /// Return the epoch the iterator instance is currently at.
    const epoch_td
    epoch() noexcept { return *m_epoch_iter;}

    /// Return the entry the iterator instance is currently at.
    const entry
    data() noexcept { return *m_data_iter; }

private:
    const timeseries<T, F>&                        m_timeseries;
    typename std::vector<entry>::const_iterator    m_data_iter;
    typename std::vector<epoch_td>::const_iterator m_epoch_iter;
}; // class timeseries_const_iterator


// @todo This class needs seriously more work (!!)
// when loping i must use hit_the_end() else it wont work!
// wtf?
template<class T, class F>
    class running_window
{
public:

    using epoch_td    = ngpt::datetime<T>;
    using interval_td = ngpt::datetime_interval<T>;
    using tflag       = ngpt::flag<F>;
    using entry       = ngpt::data_point<F>;
    using record      = std::tuple<epoch_td&, entry&>;

    explicit
    running_window(timeseries<T, F>& ts, const interval_td& w)
    : m_window{w},
      m_half_window{w},
      m_iterator_begin{ts},
      m_iterator{ts},
      m_iterator_end{ts},
      m_END{ts.end()}
    {
        assert( ts.data_pts() == ts.epochs() );
        m_half_window = split_window();
        m_END = ts.end();
    }

    void
    begin()
    {
        m_iterator_begin.begin();
        m_iterator.begin();
        m_iterator_end.begin();
        while (   (m_iterator_end != m_END)
               && (m_iterator_end.delta_time(m_iterator) < m_half_window) ) {
            ++m_iterator_end;
        }
        return;
    }

    void
    end()
    {
        m_iterator_begin.end();
        m_iterator.end();
        m_iterator_end.end();

        timeseries_iterator<T, F> it_start {m_iterator_begin.timeseries_ref()};
        it_start.begin();
        
        while (  (m_iterator_begin >= it_start)
              && (m_iterator.delta_time(m_iterator_begin) < m_half_window) ) {
            --m_iterator_begin;
        }
        return;
    }

    bool
    hit_the_end() const noexcept { return m_iterator == m_END;}
    
    // prefix
    running_window& operator++()
    {
        ++m_iterator;

        while (  (m_iterator > m_iterator_begin)
              && (m_iterator.delta_time(m_iterator_begin) > m_half_window) )
            ++m_iterator_begin;
        
        while (  (m_iterator_end != m_END)
              && (m_iterator_end.delta_time(m_iterator) < m_half_window) )
            ++m_iterator_end;
        
        return *this;
    }

    // postfix
    running_window operator++(int)
    {
        auto tmp {*this};
        this->operator++();
        return tmp;
    }

    bool
    operator==(const running_window& other) const noexcept
    { return (     m_iterator_begin == other.m_iterator_begin 
                && m_iterator == other.m_iterator
                && m_iterator_end == other.m_iterator_end ); }
    
    bool
    operator!=(const running_window& other) const noexcept
    { return !(this->operator==(other)); }

    timeseries_iterator<T, F>&
    first() noexcept
    { return m_iterator_begin; }
    
    timeseries_iterator<T, F>&
    last() noexcept
    { return m_iterator_end; }

    timeseries_iterator<T, F>
    vlast() noexcept
    { 
        auto it { m_iterator_end };
        return --it;
    }
    
    timeseries_iterator<T, F>&
    centre() noexcept
    { return m_iterator; }

    entry
    average() const noexcept
    {
        double mean{0}, sigma{0};
        int size{0};
        for (auto it = m_iterator_begin; it != m_iterator_end; ++it) {
            mean  = (it.data().value() + size*mean)/(size+1.0e0);
            sigma = (it.data().sigma() + size*mean)/(size+1.0e0);
            ++size;
        }
        return entry{mean, sigma};
    }
    
    entry
    clean_average() const noexcept
    {
        double mean{0}, sigma{0};
        int size{0};
        for (auto it = m_iterator_begin; it != m_iterator_end; ++it) {
            if ( !it.data().skip() ) {
                mean  = (it.data().value() + size*mean)/(size+1.0e0);
                sigma = (it.data().sigma() + size*mean)/(size+1.0e0);
                ++size;
            }
        }
        return entry{mean, sigma};
    }

    entry
    iqr() const noexcept
    {
        std::vector<double> vals, sigmas;
        int size = m_iterator_end.distance_from(m_iterator_begin);
        vals.reserve(size);
        sigmas.reserve(size);

        for (auto it = m_iterator_begin; it != m_iterator_end; ++it) {
            vals.push_back(it.data().value()); 
            sigmas.push_back(it.data().sigma()); 
        }

#ifdef debug
        assert( vals.size() == sigmas.size() && vals.size() == size );
#endif
        if (size == 1) return entry{vals[0], sigmas[0]};

        std::sort(vals.begin(), vals.end());
        std::sort(sigmas.begin(), sigmas.end());
        std::size_t half_size = size/2;
        if (size%2) { /* odd size */
            double vq1 = vals[half_size/2];
            double vq3 = vals[half_size+half_size/2+1];
            double sq1 = sigmas[half_size/2];
            double sq3 = sigmas[half_size+half_size/2+1];
            return entry {vq3-vq1, sq3-sq1};
        } else {
            double vq1 = vals[half_size/2];
            double vq3 = vals[half_size+half_size/2];
            double sq1 = sigmas[half_size/2];
            double sq3 = sigmas[half_size+half_size/2];
            return entry {vq3-vq1, sq3-sq1};
        }
    }

    entry
    median() const noexcept
    {
        std::vector<double> vals, sigmas;
        int size = m_iterator_end.distance_from(m_iterator_begin);
        vals.reserve(size);
        sigmas.reserve(size);

        for (auto it = m_iterator_begin; it != m_iterator_end; ++it) {
            vals.push_back(it.data().value()); 
            sigmas.push_back(it.data().sigma()); 
        }

#ifdef debug
        assert( vals.size() == sigmas.size() && vals.size() == size );
#endif

        if (size%2) { /* odd size */
            std::nth_element(vals.begin(), vals.begin() + size/2, vals.end());
            std::nth_element(sigmas.begin(), sigmas.begin() + size/2, sigmas.end());
            return entry {vals[size/2], sigmas[size/2]};
        } else {
            std::nth_element(vals.begin(), vals.begin() + size/2+1, vals.end());
            std::nth_element(sigmas.begin(), sigmas.begin() + size/2+1, sigmas.end());
            return entry {(vals[size/2]+vals[size/2-1])/2.0,
                (sigmas[size/2]+sigmas[size/2-1])/2.0};
        }
    }
    
    entry
    clean_iqr() const
    {
        std::vector<double> vals, sigmas;
        int size = m_iterator_end.distance_from(m_iterator_begin);
        vals.reserve(size);
        sigmas.reserve(size);

        for (auto it = m_iterator_begin; it != m_iterator_end; ++it) {
            if ( !it.data().skip() ) {
                vals.push_back(it.data().value()); 
                sigmas.push_back(it.data().sigma());
            }
        }
        
        size = vals.size();
        if (!size) {
            std::cerr<<"\n[ERROR] Running window with no points!";
            throw std::runtime_error("running_window::clean_iqr");
        } else if (size == 1) {
            return entry{vals[0], sigmas[0]};
        }

        std::sort(vals.begin(), vals.end());
        std::sort(sigmas.begin(), sigmas.end());
        std::size_t half_size = size/2;
        if (size%2) { /* odd size */
            double vq1 = vals[half_size/2];
            double vq3 = vals[half_size+half_size/2+1];
            double sq1 = sigmas[half_size/2];
            double sq3 = sigmas[half_size+half_size/2+1];
            return entry {vq3-vq1, sq3-sq1};
        } else {
            double vq1 = vals[half_size/2];
            double vq3 = vals[half_size+half_size/2];
            double sq1 = sigmas[half_size/2];
            double sq3 = sigmas[half_size+half_size/2];
            return entry {vq3-vq1, sq3-sq1};
        }
    }
    
    entry
    clean_median() const
    {
        std::vector<double> vals, sigmas;
        int size = m_iterator_end.distance_from(m_iterator_begin);
        vals.reserve(size);
        sigmas.reserve(size);

        for (auto it = m_iterator_begin; it != m_iterator_end; ++it) {
            if ( !it.data().skip() ) {
                vals.push_back(it.data().value()); 
                sigmas.push_back(it.data().sigma());
            }
        }

        size = vals.size();
        if (!size) {
            std::cerr<<"\n[ERROR] Running window with no points!";
            throw std::runtime_error("running_window::clean_median");
        } else if (size == 1) {
            return entry {vals[0], sigmas[0]};
        }

        if (size%2) { /* odd size */
            std::nth_element(vals.begin(), vals.begin() + size/2+1, vals.end());
            std::nth_element(sigmas.begin(), sigmas.begin() + size/2+1, sigmas.end());
            return entry {vals[size/2], sigmas[size/2]};
        } else {
            std::nth_element(vals.begin(), vals.begin() + size/2+1, vals.end());
            std::nth_element(sigmas.begin(), sigmas.begin() + size/2+1, sigmas.end());
            return entry {(vals[size/2]+vals[size/2-1])/2.0,
                (sigmas[size/2]+sigmas[size/2-1])/2.0};
        }
    }

private:
    
    /// Split the (initial, integral) window into half.
    interval_td
    split_window() noexcept
    {
        // TODO should assert that the interval is > 0
        auto mjd = m_window.days();
        modified_julian_day::underlying_type t_mjd {mjd.as_underlying_type()/2};

        auto sec = m_window.secs();
        typename T::underlying_type t_sec {sec.as_underlying_type()/2};

        if ( t_mjd % 2 ) t_sec += (T::max_in_day/2);

        interval_td half_w {modified_julian_day{t_mjd}, T{t_sec}};

        return half_w;
    }

    interval_td                m_window,
                               m_half_window;
    timeseries_iterator<T, F>  m_iterator_begin,
                               m_iterator,
                               m_iterator_end,
                               m_END;
};


template<class T, class F>
    void
    nikolaidis(timeseries<T, F>& residuals, timeseries<T, F>& original_ts,
        const datetime_interval<T>& window)
{
#ifdef DEBUG
    std::size_t num_of_outliers = 0, valid_pts = 0;
#endif
    using entry = typename timeseries<T, F>::entry;
    double res;
    entry median, iqr;
    std::size_t counter = 0;

    for (auto rw_it  = residuals.rw_begin(window);
              !rw_it.hit_the_end();
              ++rw_it)
    {
        if ( !rw_it.centre().data().skip() ) {
#ifdef DEBUG
            ++valid_pts;
#endif
            res    = rw_it.centre().data().value();
            median = rw_it.clean_median();
            iqr    = rw_it.clean_iqr();
            if ( std::abs(res-median.value()) > 3e0*iqr.value() ) {
                rw_it.centre().mark(pt_marker::outlier);
                original_ts.mark(counter, pt_marker::outlier);
#ifdef DEBUG
                ++num_of_outliers;
#endif
            }
        }
        ++counter;
    }
    assert( counter == original_ts.data_pts() );
#ifdef DEBUG
    std::cout<<"\n[DEBUG] Outliers: " << num_of_outliers << "/" 
        << valid_pts << " = ";
    printf("%5.2f", ((double)num_of_outliers/valid_pts)*100);
#endif
    return;
}

} // end namespace ngpt

#endif

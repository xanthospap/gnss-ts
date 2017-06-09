#ifndef __NGPT_TIMESERIES__
#define __NGPT_TIMESERIES__

// standard headers
#include <vector>
#include <iterator>
#include <algorithm>
#include <tuple>
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
    data_point(double val=0.0, double sigma=1.0, tflag f=tflag{})
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
/// @todo     Should the time-series be always in correct time-order. Say more
///           about this ....
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) of the time-series, will have a time-stamp of type
///           ngpt::datetime<T>.
/// @param F  An enumeration class to act as flag; each data point of the
///           time-series will have a flag of type ngpt::flag<F> (see class
///           data_point for details).
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

    /*
    std::size_t
    ts2array(F* epoch_arr, F* val_arr, F& ave, F& var) const
    {   
        std::size_t N {this->data_pts() - this->skipped_pts()};
        F prev_ave{0};
        auto ts_start {this->cbegin()},
             ts_stop  {this->cend()};
        std::size_t index {0};

        for (auto it = ts_start; it != ts_stop; ++it) {
            if ( !it.data().skip() ) {
                epoch_arr[index] = it.epoch().as_mjd();
                val_arr[index]   = it.data().value();
                prev_ave = ave;
                ave += (val_arr[index]-ave)/(index+1e0);
                var += (val_arr[index]-prev_ave)*(val_arr[index]-ave);
                ++index;
            }
        }
        var /= (N-1e0);
        return N;
    }
    */
    
    /// Constructor. If a vector of epochs is passed in, then we know we have
    /// our epochs. In this case, reserve space (memory) for the data points.
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
            if (end < start || end >= ts.m_data.size()) {
                throw std::domain_error("Invalid start/stop indexes for ts split()");
            }
            m_data.reserve(end-start);
            double sz;
            m_mean    = 0e0;
            m_skipped = 0;
            for (std::size_t i = start; i < end; ++i) {
                sz = i - start;
                m_data.emplace_back(ts[i]);
                m_mean = (ts[i].value() + sz*m_mean)/(sz+1e0);
                if ( m_data[i].skip() ) ++m_skipped;
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
    double
    add_point(entry&& e)
    {
        double sz = static_cast<double>(m_data.size());
        m_data.emplace_back(e);
        if ( e.skip() ) {
            ++m_skipped;
        } else {
            m_mean = (e.value() + sz * m_mean)/(sz+1.0);
        }
        return m_mean;
    }

    /// Add a data point; returns the new mean value.
    double
    add_point(double val, double sigma=1.0, tflag f=tflag{})
    {
        return this->add_point( entry{val, sigma, f} );
    }

    /// Mark a data point given its index.
    void
    mark(std::size_t index, F f)
    {
        tflag previous = m_data[index].flag();
        m_data[index].flag().set(f);
        if ( !__skip__(previous) && __skip__(tflag{f}) ) { ++m_skipped; }
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
        auto delta_dt = ngpt::delta_date(last_epoch(), first_epoch());
        auto central_epoch { first_epoch() };
        central_epoch += delta_dt/2;
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
        auto start_epoch = this->first_valid_epoch();
        auto delta_dt    = ngpt::delta_date(start_epoch, last_valid_epoch());
        start_epoch     += delta_dt;
        return start_epoch;
    }

    // 
    auto
    qr_ls_solve(
        ngpt::ts_model<T>&    model,
        double&               a_posteriori_std_dev,
        double sigma0         = 1e-03,
        bool   set_model_epoch= false /* set or not the mean epoch for the (input) model */
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
        x = A.colPivHouseholderQr().solve(b);
        /*
        auto sm1 = A.colPivHouseholderQr();
        x        = sm1.solve(b);
        auto R   = sm1.matrixQR().triangularView<Eigen::Upper>();
        auto S   = R.transpose();
        std::cout<<"\nCoef: "<<R(1,1);
        */
        

        // residual vector u = A*x - b; note that the reisdual vector may not
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
        model.filter_parameters();
        model.dump(std::cout);

        // Cast residuals to time-series and compute a-posteriori std. dev
        // The resulting residual ts will have the same size as the originak ts,
        // where data points marked as 'skipped' will have their original value.
        a_posteriori_std_dev = 0;
        timeseries<T, F> res {*this};
        res.epoch_ptr() = this->epoch_ptr();
        idx = counter = 0;
        double residual;
        for (std::size_t i = 0; i < epochs(); i++) {
            if ( !m_data[i].skip() ) {
                residual = u(idx)/(sigma0/m_data[i].sigma());
                data_point<F> dp { residual, m_data[i].sigma(),  m_data[i].flag() };
                a_posteriori_std_dev += residual*residual;
                res[i] = dp;
                ++idx;
            } else {
                res[i] = m_data[i];
            }
        }
        a_posteriori_std_dev = std::sqrt(a_posteriori_std_dev)
            /(double)(idx-parameters);
        std::cout<<"\nA-posteriori std. deviation: " << a_posteriori_std_dev << "(m).";

        // apply outlier detection algorithm and mark them
        datetime_interval<T> window {modified_julian_day{90}, T{0}};
        nikolaidis(res, *this, window);

        return res;
    }

    /// Return a const iterator to the first entry of the data points vector
    typename std::vector<entry>::const_iterator
    citer_start() const noexcept
    { return m_data.cbegin(); }

    /// Return an (non-const) iterator to the first entry of the data points vector
    typename std::vector<entry>::iterator
    iter_start() noexcept
    { return m_data.begin(); }
    
    /// Return a const iterator to the last+1 (=end) entry of the data points vector
    typename std::vector<entry>::const_iterator
    citer_stop() const noexcept
    { return m_data.cend(); }
    
    /// Return a (non-const) iterator to the last+1 (=end) entry of the data points vector
    typename std::vector<entry>::iterator
    iter_stop() noexcept
    { return m_data.end(); }

    timeseries_iterator<T, F>
    begin()
    {
        return timeseries_iterator<T, F> {*this};
    }

    timeseries_iterator<T, F>
    end()
    {
        timeseries_iterator<T, F> it {*this};
        it.end();
        return it;
    }
    
    timeseries_const_iterator<T, F>
    cbegin()
    const
    {
        return timeseries_const_iterator<T, F> {*this};
    }

    timeseries_const_iterator<T, F>
    cend()
    const
    {
        timeseries_const_iterator<T, F> it {*this};
        it.end();
        return it;
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

    explicit
    timeseries_iterator(timeseries<T, F, typename std::enable_if_t<T::is_of_sec_type>>& ts)
    : m_timeseries{ts},
      m_data_iter{ts.iter_start()},
      m_epoch_iter{ts.epoch_ptr()->begin()}
    {
        assert( ts.epoch_ptr() && ts.data_pts() == ts.epochs() );
    }

    //
    timeseries_iterator(const timeseries_iterator& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {}

    timeseries_iterator(timeseries_iterator&& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {
    }

    //
    timeseries_iterator& operator=(const timeseries_iterator& it) noexcept
    {
        if (*this != it) {
            m_timeseries = it.m_timeseries;
            m_data_iter  = it.m_data_iter;
            m_epoch_iter = it.m_epoch_iter;
        }
        return *this;
    }
    
    timeseries_iterator& operator=(timeseries_iterator&& it) noexcept
    {
        if (*this != it) {
            m_timeseries = std::move(it.m_timeseries);
            m_data_iter  = std::move(it.m_data_iter);
            m_epoch_iter = std::move(it.m_epoch_iter);
        }
        return *this;
    }

    void begin() noexcept
    {
        m_data_iter  = m_timeseries.iter_start();
        m_epoch_iter = m_timeseries.epoch_ptr()->begin();
    }

    void end() noexcept
    {
        m_data_iter  = m_timeseries.iter_stop();
        m_epoch_iter = m_timeseries.epoch_ptr()->end();
    }

    void advance() noexcept
    {
        ++m_data_iter;
        ++m_epoch_iter;
    }

    std::size_t
    index() noexcept
    {
        return std::distance(m_timeseries.iter_start(), m_data_iter);
    }

    bool
    operator==(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter == it.m_data_iter && m_epoch_iter == it.m_epoch_iter );
    }

    bool
    operator!=(const timeseries_iterator& it) const noexcept
    {
        return !((*this) == it);
    }

    bool
    operator>(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter > it.m_data_iter && m_epoch_iter > it.m_epoch_iter );
    }
    
    bool
    operator>=(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter >= it.m_data_iter && m_epoch_iter >= it.m_epoch_iter );
    }
    
    bool
    operator<(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter < it.m_data_iter && m_epoch_iter < it.m_epoch_iter );
    }
    
    bool
    operator<=(const timeseries_iterator& it) const noexcept
    {
        return ( m_data_iter <= it.m_data_iter && m_epoch_iter <= it.m_epoch_iter );
    }

    timeseries_iterator /* pointer arithmetic */
    operator+(int n) const noexcept
    {
        timeseries_iterator ti {*this};
        ti.m_data_iter + n;
        ti.m_epoch_iter + n;
        return ti;
    }
    timeseries_iterator /* pointer arithmetic */
    operator-(int n) const noexcept
    {
        timeseries_iterator ti {*this};
        ti.m_data_iter - n;
        ti.m_epoch_iter - n;
        return ti;
    }

    // prefix
    timeseries_iterator& operator++()
    {
        ++m_data_iter;
        ++m_epoch_iter;
        return *this;
    }
    
    // prefix
    timeseries_iterator& operator--()
    {
        --m_data_iter;
        --m_epoch_iter;
        return *this;
    }

    // postfix
    timeseries_iterator operator++(int)
    {
        auto tmp {*this};
        this->operator++();
        return tmp;
    }

    int
    distance_from(const timeseries_iterator& it) const
    {
        assert( &m_timeseries == &it.m_timeseries );
        return std::distance(it.m_data_iter, this->m_data_iter);
    }

    interval_td
    delta_time(const timeseries_iterator& it) const
    {
        assert( &m_timeseries == &(it.m_timeseries) );
        return m_epoch_iter->delta_date( *(it.m_epoch_iter) );
    }

    epoch_td&
    epoch() noexcept { return *m_epoch_iter;}

    entry&
    data() noexcept { return *m_data_iter; }
    
    void
    mark(F f)
    {
        tflag previous = m_data_iter->flag();
        m_data_iter->flag().set(f);
        if ( !__skip__(previous) && __skip__(tflag{f}) ) {
            m_timeseries.skipped_pts() = m_timeseries.skipped_pts() + 1;
        }
    }

    timeseries<T, F>&
    timeseries_ref() noexcept
    { return m_timeseries; }

private:
    timeseries<T, F>&                        m_timeseries;
    typename std::vector<entry>::iterator    m_data_iter;
    typename std::vector<epoch_td>::iterator m_epoch_iter;
};

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

    explicit
    timeseries_const_iterator(const timeseries<T, F, typename std::enable_if_t<T::is_of_sec_type>>& ts)
    : m_timeseries{ts},
      m_data_iter{ts.citer_start()},
      m_epoch_iter{ts.epoch_ptr()->cbegin()}
    {
        assert( ts.epoch_ptr() && ts.data_pts() == ts.epochs() );
    }

    //
    timeseries_const_iterator(const timeseries_const_iterator& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {}

    timeseries_const_iterator(timeseries_const_iterator&& it) noexcept
    : m_timeseries{it.m_timeseries},
      m_data_iter{it.m_data_iter},
      m_epoch_iter{it.m_epoch_iter}
    {
    }

    //
    timeseries_const_iterator& operator=(const timeseries_const_iterator& it) noexcept
    {
        if (*this != it) {
            m_timeseries = it.m_timeseries;
            m_data_iter  = it.m_data_iter;
            m_epoch_iter = it.m_epoch_iter;
        }
        return *this;
    }
    
    timeseries_const_iterator& operator=(timeseries_const_iterator&& it) noexcept
    {
        if (*this != it) {
            m_timeseries = std::move(it.m_timeseries);
            m_data_iter  = std::move(it.m_data_iter);
            m_epoch_iter = std::move(it.m_epoch_iter);
        }
        return *this;
    }

    void begin() noexcept
    {
        m_data_iter  = m_timeseries.citer_start();
        m_epoch_iter = m_timeseries.epoch_ptr()->cbegin();
    }

    void end() noexcept
    {
        m_data_iter  = m_timeseries.citer_stop();
        m_epoch_iter = m_timeseries.epoch_ptr()->cend();
    }

    void advance() noexcept
    {
        ++m_data_iter;
        ++m_epoch_iter;
    }

    std::size_t
    index() noexcept
    {
        return std::distance(m_timeseries.citer_start(), m_data_iter);
    }

    bool
    operator==(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter == it.m_data_iter && m_epoch_iter == it.m_epoch_iter );
    }

    bool
    operator!=(const timeseries_const_iterator& it) const noexcept
    {
        return !((*this) == it);
    }

    bool
    operator>(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter > it.m_data_iter && m_epoch_iter > it.m_epoch_iter );
    }
    
    bool
    operator>=(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter >= it.m_data_iter && m_epoch_iter >= it.m_epoch_iter );
    }
    
    bool
    operator<(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter < it.m_data_iter && m_epoch_iter < it.m_epoch_iter );
    }
    
    bool
    operator<=(const timeseries_const_iterator& it) const noexcept
    {
        return ( m_data_iter <= it.m_data_iter && m_epoch_iter <= it.m_epoch_iter );
    }

    timeseries_const_iterator /* pointer arithmetic */
    operator+(int n) const noexcept
    {
        timeseries_const_iterator ti {*this};
        ti.m_data_iter  += n;
        ti.m_epoch_iter += n;
        return ti;
    }

    timeseries_const_iterator /* pointer arithmetic */
    operator-(int n) const noexcept
    {
        timeseries_const_iterator ti {*this};
        ti.m_data_iter  -= n;
        ti.m_epoch_iter -= n;
        return ti;
    }

    // prefix
    timeseries_const_iterator& operator++()
    {
        ++m_data_iter;
        ++m_epoch_iter;
        return *this;
    }
    
    // prefix
    timeseries_const_iterator& operator--()
    {
        --m_data_iter;
        --m_epoch_iter;
        return *this;
    }

    // postfix
    timeseries_const_iterator operator++(int)
    {
        auto tmp {*this};
        this->operator++();
        return tmp;
    }

    int
    distance_from(const timeseries_const_iterator& it) const
    {
        assert( &m_timeseries == &it.m_timeseries );
        return std::distance(it.m_data_iter, this->m_data_iter);
    }

    interval_td
    delta_time(const timeseries_const_iterator& it) const
    {
        assert( &m_timeseries == &(it.m_timeseries) );
        return m_epoch_iter->delta_date( *(it.m_epoch_iter) );
    }

    epoch_td
    epoch() noexcept { return *m_epoch_iter;}

    entry
    data() noexcept { return *m_data_iter; }

private:
    const timeseries<T, F>&                        m_timeseries;
    typename std::vector<entry>::const_iterator    m_data_iter;
    typename std::vector<epoch_td>::const_iterator m_epoch_iter;
};


// TODO
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
    clean_iqr() const noexcept
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
        if (size == 1) {
            std::cout<<"------ZERO SIZE VECTOR--------------";
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
    clean_median() const noexcept
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
        if (size<5) std::cout<<"\nToo few points:";

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
    nikolaidis(timeseries<T, F>& residuals, timeseries<T, F>& original_ts, const datetime_interval<T>& window)
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
            // auto foo = rw_it.centre().epoch();
            // std::cout<<"\n"<<foo.as_mjd()<<" " <<original_ts[counter].value()<<" "<<res;
            if ( std::abs(res-median.value()) > 3.0*iqr.value() ) {
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
    std::cout<<"\nOutliers: " << num_of_outliers << "/" << valid_pts << " = " << ((double)num_of_outliers/valid_pts)*100;
#endif
    return;
}

/*
template<class T, class F>
    void
    three_sigma(timeseries<T, F>& residuals, timeseries<T, F>& original_ts, double sigma)
{
#ifdef DEBUG
    std::size_t num_of_outliers = 0, valid_pts = 0;
#endif
    double res;
    std::size_t counter = 0;

    for (auto rw_it  = residuals.begin();
              rw_it != residuals.end();
              ++rw_it)
    {
        if ( !rw_it.data().skip() ) {
#ifdef DEBUG
            ++valid_pts;
#endif
            res    = rw_it.data().value();
            //median = rw_it.clean_median();
            //iqr    = rw_it.clean_iqr();
            if ( std::abs(res) > 3.0*sigma ) {
                rw_it.data().flag().set(pt_marker::outlier);
                // std::cout<<"\nmarking cause: |" << res << "-"<< median.value() << "| > 3.0 * " << iqr.value();
                original_ts[counter].flag().set(pt_marker::outlier);
#ifdef DEBUG
                ++num_of_outliers;
#endif
            }
        }
        ++counter;
    }
    assert( counter == original_ts.data_pts() );
#ifdef DEBUG
    std::cout<<"\nOutliers: " << num_of_outliers << "/" << valid_pts << " = " << ((double)num_of_outliers/valid_pts)*100;
#endif
    return;
}*/

} // end namespace ngpt

#endif

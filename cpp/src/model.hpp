#ifndef __NGPT_TS_MODEL__
#define __NGPT_TS_MODEL__

#include <iostream>
#ifdef DEBUG
#include <cfenv>
#include <cstdlib>
#endif
// Eigen headers
#include "eigen3/Eigen/Core"
// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"
// gtms headers
#include "psd.hpp"
#include "event_list.hpp"

namespace ngpt
{

/// A class to represent a "jump" (i.e. offset) of a time-series. A 'jump' is
/// determined by the value of the offset and the epoch (time tag) at which it
/// occured.
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) will have a time-stamp of type ngpt::datetime<T>.
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_jump
{
public:
    /// Constructor.
    /// @param[in] start  The epoch the jump happened.
    /// @param[in] offset The value of the offset (default is 0).
    explicit
    md_jump(ngpt::datetime<T> start, double offset=0e0, double stdd=0e0) noexcept
    : m_start {start},
      m_offset{offset},
      m_stddev{stdd}
    {};

    /// Get the epoch the jump happened at (const).
    datetime<T>
    start() const noexcept { return m_start; }

    /// Get the value of the offset (const).
    double
    value() const noexcept { return m_offset; }
    
    /// Get the value of the offset (non-const).
    double&
    value() noexcept { return m_offset; }
    
    /// Get the std. deviation value of the offset (const).
    double
    stddev() const noexcept { return m_stddev; }
    
    /// Get the value of the offset (non-const).
    double&
    stddev() noexcept { return m_stddev; }

private:
    ngpt::datetime<T> m_start;  ///< When the "jump" occured.
    double            m_offset; ///< Value (i.e. offset amplitude).
    double            m_stddev; ///< Std. deviation of the estimated value.

}; // md_jump

/// A class to represent a harmonic signal (in a time-series). The harmonic
/// signal is described by in and out-of-phase amplitudes, a period/frequency,
/// and the starting and ending times (i.e. its validity interval).
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) will have a time-stamp of type ngpt::datetime<T>.
///@note      Angular frequency (i.e. ω) is given by: ω = 2πf = ω = 2π/T
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_harmonics
{
public:
    /// Constructor.
    /// @param[in] period    The period of the harmonic.
    /// @param[in] start     Start of the validity interval (default value
    ///                      datetime<T>::min()).
    /// @param[in] stop      End of the validity interval (default value
    ///                      datetime<T>::max()).
    /// @param[in] in_phase  Amplitude of the in-phase component.
    /// @param[in] out_phase Amplitude of the out-of-phase component.
    //
    explicit
    md_harmonics(double period,
        ngpt::datetime<T> start=ngpt::datetime<T>::min(),
        ngpt::datetime<T> stop = ngpt::datetime<T>::max(),
        double in_phase=0e0, double out_phase=0e0,
        double in_phase_stddev=0e0, double out_phase_stddev=0e0) noexcept
    : m_start    {start},
      m_stop     {stop},
      m_afreq    {D2PI/period},
      m_in_phase {in_phase},
      m_in_stddev{in_phase_stddev},
      m_out_phase{out_phase},
      m_out_stddev{out_phase_stddev}
    {};

    /// Get the starting epoch (of the validity interval).
    ngpt::datetime<T>
    start() const noexcept { return m_start; }

    /// Get the ending epoch (of the validity interval).
    ngpt::datetime<T>
    stop() const noexcept { return m_stop; }

    /// Get the angular frequency (ω = 2πf).
    double
    angular_frequency() const noexcept { return m_afreq; }

    /// Get the period (i.e. T = 2π/ω).
    double
    period() const noexcept { return D2PI/m_afreq;}

    /// Get the in-phase component amplitude (const version).
    double
    in_phase() const noexcept { return m_in_phase; }
    
    /// Get the out-of-phase component amplitude (const version).
    double
    out_of_phase() const noexcept { return m_out_phase; }
    
    /// Get the in-phase component amplitude (non-const version).
    double&
    in_phase() noexcept { return m_in_phase; }
    
    /// Get the out-of-phase component amplitude (non-const version).
    double&
    out_of_phase() noexcept { return m_out_phase; }
    
    /// Get the in-phase component std. deviation (const version).
    double
    in_phase_stddev() const noexcept { return m_in_stddev; }
    
    /// Get the out-of-phase component std. deviation (const version).
    double
    out_phase_stddev() const noexcept { return m_out_stddev; }
    
    /// Get the in-phase component std. deviation (non-const version).
    double&
    in_phase_stddev() noexcept { return m_in_stddev; }
    
    /// Get the out-of-phase component std. deviation (non-const version).
    double&
    out_phase_stddev() noexcept { return m_out_stddev; }

    /// Get the amplitude of the harmonic signal, i.e. if the harmonic is:
    /// A*sin(ωt) + B*cos(ωt), then its aplitude is sqrt(A^2 + B^2).
    double
    amplitude() const noexcept
    { return std::sqrt(m_in_phase*m_in_phase+m_out_phase*m_out_phase); }

private:
    ngpt::datetime<T> m_start,      ///< Starting epoch of the harmonic
                      m_stop;       ///< Ending epoch of the harmonic
    double            m_afreq,      ///< Angular frequency (i.e. ω)
                      m_in_phase,   ///< In-Phase component (amplitude)
                      m_in_stddev,  ///< Std. deviation of In-Phase component
                      m_out_phase,  ///< Out-Of-Phase component (amplitude)
                      m_out_stddev; ///< Std. deviation of Out-Of-Phase component

}; // harmonic_coef

/// A class to represent a velocity change (in a time-series). The class is
/// dead simple, just holding a start and stop date (i.e. the validity interval
/// of this velocity change) and the magnitude of the velocity (change).
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) will have a time-stamp of type ngpt::datetime<T>.
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_velocity_change
{
public:

    /// Constructor. The stop date defaults to ngpt::datetime<T>::max().
    /// @param[in] start   The starting epoch of this velocity change.
    /// @param[in] stop    The end epoch of this velocity change; default is
    ///                    ngpt::datetime<T>::max(), i.e. this velocity change
    ///                    starts at epoch start and never stops.
    /// @param[in] new_vel The magnitude of the velocity change.
    explicit
    md_velocity_change(ngpt::datetime<T> start,
        ngpt::datetime<T> stop = ngpt::datetime<T>::max(),
        double new_vel=0e0, double new_vel_stddev=0e0) noexcept
    : m_start {start},
      m_stop  {stop},
      m_newvel{new_vel},
      m_newvel_stddev{new_vel_stddev}
    {};

    /// Get the start date (start of validity interval), const version.
    ngpt::datetime<T>
    start() const noexcept { return m_start; }

    /// Get the stop date (end of validity interval), const version.
    ngpt::datetime<T>
    stop() const noexcept { return m_stop; }

    /// Get the velocity change magnitude (const version).
    double
    value() const noexcept { return m_newvel; }
    
    /// Get the velocity change magnitude (non-const version).
    double&
    value() noexcept { return m_newvel; }
    
    /// Get the velocity change std. deviation (const version).
    double
    stddev() const noexcept { return m_newvel_stddev; }
    
    /// Get the velocity change std. deviation (non-const version).
    double&
    stddev() noexcept { return m_newvel_stddev; }
private:

    ngpt::datetime<T> m_start,  ///< Start of validity niterval.
                      m_stop;   ///< End of validity interval.
    double            m_newvel, ///< Magnitude of velocity change.
                      m_newvel_stddev; ///< Velocity change std. deviation.

}; // md_velocity_change

/// This class represents a (time-series) model.
///
/// The purpose of this class is to describe a (mathematical) model of a
/// time-series. So it should serve several requirements, e.g.
/// - given an epoch (i.e. datetime), compute the value of the model,
/// - be used within Least-Squares modules to perform adjustments (i.e. fit)
/// The model can consist of any number of md_jump, md_harmonics and
/// md_velocity_change, plus a simple linear model (i.e. x0 and Vx0 parameters);
/// so if no md_jump, md_harmonics or md_velocity_change are added, then the
/// model defaults to a linear one.
/// Every instance also hold am epoch (i.e. datetime<T>) member variable,
/// which is the central epoch of computation/estimation.
template<class T>
    class ts_model
{
public:

    /// Constructor; this default to a linear model.
    ts_model() noexcept : m_x0{0e0}, m_vx{0e0} {};

    /// A model is non-linear if it must estimate PSD/earthquakes.
    bool
    is_linear() const noexcept { 
        if ( m_earthqs.size() ) {
            for (auto j = m_earthqs.cbegin(); j != m_earthqs.cend(); ++j)
                if ( j->parameters() != 1 )
                    return false;
        }
        return true;
    }
    
    /// Constructor using an event_list instance.
    /// @param[in] events An instance of type event_list; all events recorded
    ///                   in this instance, are going to be included in the
    ///                   constructed model. This includes:
    ///                   - jumps
    ///                   - velocity changes
    ///                   - earthquakes
    explicit
    ts_model(const event_list<T>& events) noexcept
    : m_x0{0e0},
      m_vx{0e0}
    {
      for (auto it = events.it_cbegin(); it != events.it_cend(); ++it)
      {
        if ( it->event_type() == ts_event::jump ) {
          m_jumps.emplace_back(it->epoch());
        } else if ( it->event_type() == ts_event::velocity_change ) {
          m_vel_changes.emplace_back(it->epoch());
        } else if ( it->event_type() == ts_event::earthquake ) {
          m_earthqs.emplace_back(it->epoch(), psd_model::pwl);
        }
      }
      ngpt::datetime<T> t0 {modified_julian_day{0}, T{0}};
      m_mean_epoch = t0;
    }

    /// Given a vector of epochs (aka datetime<T> instances), compute the model
    /// values at the corresponding datetimes.
    /// @param[in] epochs A (reference to) vector of epochs (aka datetime<T>)
    /// @return    A vector of doubles, where each entry i is the corresponding
    ///            model value at epochs[i].
    /// @note      
    ///          - The central epoch for the computation is the instance's
    ///            mean epoch (aka m_mean_epoch).
    ///          - Delta-times (Δt) are considered to the differences of each
    ///            epoch from the mean epoch as mjd (aka difference is in days)
    /// @todo    I do not like the fact that the velocity model parameters are
    ///          set to annual.
    std::vector<double>
    make_model(const std::vector<datetime<T>>& epochs)
    const
    {
        assert( epochs.size() );
        double value    = 0e0,
               mean_mjd = m_mean_epoch.as_mjd(),
               dt       = 0e0;
        std::vector<double> vals;
        vals.reserve( epochs.size() );

        for (const auto& t : epochs) {

            // compute δt in days
            dt    = t.as_mjd() - mean_mjd;
            
            // coef for constant (linear) velocity i.e. m/year
            value = m_x0 + m_vx*(dt/365.25);

            // Harmonic coefficients for each period ...
            for (auto j = m_harmonics.cbegin(); j != m_harmonics.cend(); ++j) {
                if (t >= j->start() && t < j->stop() ) {
                    // cosinus or in-phase
                    value += j->in_phase()*std::cos(j->angular_frequency()*dt);
                    // sinus or out-of-phase
                    value += j->out_of_phase()*std::sin(j->angular_frequency()*dt);
                }
            }

            // Set up jumps ...
            for (auto j = m_jumps.cbegin(); j != m_jumps.cend(); ++j) {
                if ( t >= j->start() )
                    value += j->value();
            }

            // Set up velocity changes ...
            for (auto j = m_vel_changes.cbegin(); j != m_vel_changes.cend(); ++j) {
                if (t >= j->start() && t < j->stop()) {
                    auto   dti  = ngpt::delta_date(t, j->start());
                    double dtj  = dti.as_mjd();
                    value += j->value()*(dtj/365.25);
                }
            }
            
            // Set up earthquake PSD corrections
            for (auto j = m_earthqs.cbegin(); j != m_earthqs.cend(); ++j) {
                if ( t >= j->start() ) {
                    /*
                    auto   dti = ngpt::delta_date(t, j->start());
                    double dtj = dti.as_mjd()/365.25e0;
                    double t1  = (dtj)/j->tau_log();
                    double t2  = (dtj)/j->tau_exp();
                    value += j->mag_log() * std::log(1e0+t1);
                    value += j->mag_exp() * (1e0-std::exp(-t2));
                    */
                    value += j->value_at(t);
/*
#ifdef DEBUG
                    if ( std::fetestexcept(FE_ALL_EXCEPT & ~FE_INEXACT) ) {
                        std::cerr<<"\n\nFloating Point Exception at \"make_model()\"";
                        std::cerr<<"\nHere are the variables:";
                        std::cerr<<"\n\tdti.as_mjd()="<<dti.as_mjd();
                        std::cerr<<"\n\ttau_log     ="<<j->tau_log();
                        std::cerr<<"\n\ttau_exp     ="<<j->tau_exp();
                        std::cerr<<"\n\tdtj         ="<<dtj;
                        std::cerr<<"\n\tt1          ="<<t1;
                        std::cerr<<"\n\tt2          ="<<t2;
                        std::cerr<<"\n\tmag_log     ="<<j->mag_log();
                        std::cerr<<"\n\tmag_exp     ="<<j->mag_exp();
                        std::cerr<<"\n\tt.as_mjd()  ="<<t.as_mjd();
                        std::cerr<<"\n\tindex       ="<<vals.size()+1;
                        std::cerr<<"\nFP Exceptions Set:";
                        std::cerr<<"\n\tFE_DIVBYZERO: "<<std::fetestexcept(FE_DIVBYZERO);
                        std::cerr<<"\n\tFE_INEXACT  : "<<std::fetestexcept(FE_INEXACT);
                        std::cerr<<"\n\tFE_INVALID  : "<<std::fetestexcept(FE_INVALID);
                        std::cerr<<"\n\tFE_OVERFLOW : "<<std::fetestexcept(FE_OVERFLOW);
                        std::cerr<<"\n\tFE_UNDERFLOW: "<<std::fetestexcept(FE_UNDERFLOW);
                        std::exit(1);
                    }
#endif
*/
                }
            }

            vals.emplace_back(value);
        }
        return vals;
    }

    /// Given a start an end date and a delta-time (aka a datetime grid),
    /// compute the model value for each of the datetime points, i.e. for
    /// each epoch in the interval [start, stop) with step = tinterval.
    /// @param[in] start     The starting epoch, first epoch to get the model
    ///                      values for.
    /// @param[in] stop      The ending epoch; note that the ending value is not
    ///                      in the range, i.e. the range is [start, stop)
    /// @param[in] tinterval The step size to construct epochs in the range
    ///                      [start, stop)
    /// @param[in] epoch_vec Pointer to a vector of epochs (aka datetime<T>).
    ///                      If not provided (defaults to nullptr), then it is
    ///                      not used; else, it is cleared and filled in with
    ///                      the time points for which the model value was
    ///                      computed for. I.e. it will contain all points
    ///                      within the interval [start, stop) with
    ///                      step = tinterval.
    /// @return              A vector of doubles where each entry is the
    ///                      corresponding model value for the interval
    ///                      [start, stop). I.e. result[i] is the model value
    ///                      for time t = start + i*tinterval.
    std::vector<double>
    make_model(datetime<T> start, datetime<T> stop,
        datetime_interval<T> tinterval,
        std::vector<datetime<T>>* epoch_vec=nullptr)
    const
    {
        assert( start < stop );
        assert( tinterval.as_mjd() > 0e0 );
        
        // get size hint
        std::size_t size_hint = (std::size_t)((stop.as_mjd() - start.as_mjd())
                                / tinterval.as_mjd()) + 1;

        // construct vector of epochs or fill in the one passed in
        std::vector<datetime<T>>* tvec;
        std::vector<datetime<T>> epochs;
        if ( epoch_vec ) {
            epoch_vec->clear();
            if ( epoch_vec->capacity() < size_hint )
                epoch_vec->reserve(size_hint);
            tvec = epoch_vec;
        } else {
            epochs.reserve(size_hint);
            tvec = &epochs;
        }
        datetime<T> t = start;
        while ( t < stop ) {
            tvec->emplace_back( t );
            t += tinterval;
        }

        // let another function do the dirty work.
        return this->make_model(*tvec);
    }
    
    /// Given an Eigen vector, this function will translate all values of the
    /// vector to the corresponding model parameters. The size of the vector
    /// (i.e. num of rows) must match the number of parameters of the instance.
    /// This function is used e.g. in the following scenario:
    /// - construct a model instance (to fit to some time-series)
    /// - use the model to fit the time-series (e.g. via least squares)
    /// - assign the estimated parameters to the model parameters via this
    ///   function
    /// @param[in] x_estim An Eigen vector with the values of the model
    ///                    parameters
    /// @note If the model is non-linear, then the estimates (from the input
    ///       solution vector) are going to be added to the model parameters
    ///       (not assigned).
    void
    assign_solution_vector(const Eigen::VectorXd& x_estim, const Eigen::MatrixXd*Q = nullptr)
    {
        if ( !this->is_linear() ) {
            return this->add_solution_vector(x_estim);
        }

        assert( x_estim.size() >= 2 
                && (int)x_estim.size() == (int)this->parameters() );

        std::size_t idx = 0;
        m_x0 = x_estim(idx); ++idx;
        m_vx = x_estim(idx); ++idx;
        for (auto j = m_harmonics.begin(); j != m_harmonics.end(); ++j) {
            j->in_phase()     = x_estim(idx); ++idx;
            j->out_of_phase() = x_estim(idx); ++idx;
        }
        for (auto j = m_jumps.begin(); j != m_jumps.end(); ++j) {
            j->value() = x_estim(idx); ++idx;
        }
        for (auto j = m_vel_changes.begin(); j != m_vel_changes.end(); ++j) {
            j->value() = x_estim(idx); ++idx;
        }
        // in the linear case, PSDs are only one-parameter models
        for (auto j = m_earthqs.begin(); j != m_earthqs.end(); ++j) {
            j->a1() = x_estim(idx); ++idx;
        }

        if ( Q ) {
            assert( Q->rows() == Q->cols() && Q->rows() == (int)this->parameters() );
            idx = 0;
            m_x0_stddev = (*Q)(idx, idx); ++idx;
            m_vx_stddev = (*Q)(idx, idx); ++idx;
            for (auto j = m_harmonics.begin(); j != m_harmonics.end(); ++j) {
                j->in_phase_stddev()  = (*Q)(idx, idx); ++idx;
                j->out_phase_stddev() = (*Q)(idx, idx); ++idx;
            }
            for (auto j = m_jumps.begin(); j != m_jumps.end(); ++j) {
                j->stddev() = (*Q)(idx, idx); ++idx;
            }
            for (auto j = m_vel_changes.begin(); j != m_vel_changes.end(); ++j) {
                j->stddev() = (*Q)(idx, idx); ++idx;
            }
            // in the linear case, PSDs are only one-parameter models
            for (auto j = m_earthqs.begin(); j != m_earthqs.end(); ++j) {
                j->a1() = x_estim(idx); ++idx;
            }
        }

        return;
    }

    /// Write an instance to an output stream.
    /// @todo 
    ///         - There should be a corresponding read function
    ///         - Why is this not the '<<' operator?
    void
    dump(std::ostream& os) const
    {
        // write mean epoch (of computation)
        if (m_mean_epoch.mjd() != modified_julian_day{0})
            os << "Central Epoch      : " << strftime_ymd_hms<T>(m_mean_epoch)<<"\n";

        // write linear coefficients
        os << "Linear Coefficients:";
        os << "\n\t" << m_x0 << "\n\t" << m_vx;
        
        // write harmonics (if any)
        if ( m_harmonics.size() ) os << "\nHarmonic Coefficients:";
        for (auto j = m_harmonics.begin(); j != m_harmonics.end(); ++j) {
            os << "\n\t" << j->in_phase()     << " (in-phase)   ";
            os << "From: " << (j->start()==datetime<T>::min()?"start":strftime_ymd_hms<T>(j->start()) );
            os << " To: "  << (j->stop()==datetime<T>::max()?"end":strftime_ymd_hms<T>(j->stop()) );
            os << " Period: " << j->period();
            os << "\n\t" << j->out_of_phase() << " (out-of-pase)";
        }
        
        // write jumps (if any)
        if ( m_jumps.size() ) os << "\nJumps/Offsets:";
        for (auto j = m_jumps.begin(); j != m_jumps.end(); ++j) {
            os << "\n\t" << j->value() << " at " << strftime_ymd_hms<T>(j->start());
        }

        // write velocity changes (if any)
        if ( m_vel_changes.size() ) os << "\nVelocity Changes:";
        for (auto j = m_vel_changes.begin(); j != m_vel_changes.end(); ++j) {
            os << "\n\t" << j->value();
            os << " From: " << (j->start()==datetime<T>::min()?"start":strftime_ymd_hms<T>(j->start()) );
            os << " To: "  << (j->stop()==datetime<T>::max()?"end":strftime_ymd_hms<T>(j->stop()) );
        }
        
        // write earthquakes (if any)
        if ( m_earthqs.size() ) os << "\nPSD (Earthquake Post Deformation):";
        for (auto j = m_earthqs.begin(); j != m_earthqs.end(); ++j) {
            std::size_t p = j->parameters();
            os << "\nTime of earthquake: "<<strftime_ymd_hms<T>(j->start())<<" (MJD "<<j->start().as_mjd()<<")";
            if ( p == 1 ) {
                os << "\n\tOffset: "<< j->a1();
            } 
            if ( p > 1 ) {
                os << "\n\tMagnitude: "<<j->a1()<<" Tau: "<<j->t1()<<" Mdl: "<<_psd2int_(j->psd_type());
            }
            if ( p > 3 ) {
                os << "\n\tMagnitude: "<<j->a2()<<" Tau: "<<j->t2()<<" Mdl: "<<_psd2int_(j->psd_type());
            }
        }

        // all done
        return;
    }

    /// Fileter out parameters (i.e. erase them from this instance) given a
    /// minimum (absolute value). The parametrs to be considered here, are:
    /// - jumps (i.e. offsets) and
    /// - harmonics
    /// @param[in] min_abs_jump Minimum accepted absolute value for any 'jump'
    ///                         parameter. If the model contains any jump with
    ///                         amplitude less than min_abs_jump (in absolute
    ///                         value, then this jump (i.e. the md_jump instance
    ///                         is removed from the model).
    /// @param[in] min_abs_harm Minimum accepted absolute value for any 'harmonic'
    ///                         parameter. If the model contains any harmonic with
    ///                         amplitude less than min_abs_harm (in absolute
    ///                         value, then this signal (i.e. the md_harmonics
    ///                         instance is removed from the model). If the
    ///                         harmonic has the form:
    ///                         A*sin(ωt) + B*cos(ωt), then its aplitude is
    ///                         sqrt(A^2 + B^2).
    void
    filter_parameters(double min_abs_jump=1e-3, double min_abs_harm=1e-3)
    {
        typename std::vector<md_jump<T>>::iterator jit;
        while ( (jit = std::find_if(m_jumps.begin(), m_jumps.end(),
                        [=](const md_jump<T>& j){ 
                            return std::abs(j.value())<min_abs_jump;}))
                 != m_jumps.end() )
        {
            m_jumps.erase( jit );
        }
        
        typename std::vector<md_harmonics<T>>::iterator hit;
        while ( (hit = std::find_if(m_harmonics.begin(), m_harmonics.end(),
                        [=](const md_harmonics<T>& j){
                            return j.amplitude()<min_abs_harm;}))
                != m_harmonics.end() )
        {
            m_harmonics.erase( hit );
        }

        return;
    }

    /// @todo this needs documentation
    void
    dump_as_json(std::ostream& os) const
    {
        double min_mjd = datetime<T>::min().as_mjd();
        double max_mjd = datetime<T>::max().as_mjd();

        os << "{\n";

        os << "\"reference_epoch\":";
        if (m_mean_epoch.mjd() != modified_julian_day{0})
            os << m_mean_epoch.as_mjd();
        else
            os << "null";

        os << ",\n\"const_term\":" << m_x0 
           << ",\n\"velocity\":" << m_vx;
        
        if ( m_harmonics.size() ) {
            os << ",\n\"harmonics\": [";
            for (auto j = m_harmonics.begin(); j != m_harmonics.end(); ++j) {
                os << "\n{\"in_phase\":" << j->in_phase();
                os << ",\n\"from\":" << (j->start()==datetime<T>::min()?min_mjd:(j->start().as_mjd()));
                os << ",\n\"to\":"  << (j->stop()==datetime<T>::max()?max_mjd:(j->stop().as_mjd()));
                os << ",\n\"period\":" << j->period();
                os << ",\n\"out_of_phase\":" << j->out_of_phase()<<"}";
                if (j != m_harmonics.end()-1) os << ",";
            }
            os << "\n]";
        }
        
        if ( m_jumps.size() ) {
            os << ",\n\"jumps\": [";
            for (auto j = m_jumps.begin(); j != m_jumps.end(); ++j) {
                os << "\n{\"value\": "<<j->value()<<", \"at\":"<<j->start().as_mjd()<<"}";
                if (j != m_jumps.end()-1) os << ",";
            }
            os << "\n]";
        }

        if ( m_vel_changes.size() ) {
            os << ",\n\"velocity_changes\": [";
            for (auto j = m_vel_changes.begin(); j != m_vel_changes.end(); ++j) {
                os << "\n{\"value\":" << j->value() << ",";
                os << ",\n\"from\":" << (j->start()==datetime<T>::min()?min_mjd:j->start().as_mjd());
                os << ",\n\"to\":" << (j->stop()==datetime<T>::max()?max_mjd:j->stop().as_mjd())<<"}";
                if (j != m_vel_changes.end()-1) os << ",";
            }
            os << "\n]";
        }

        os << "}"; // close object

        return;
    }
    
    /// This function adds new md_harmonics instances to this (model)
    /// instance. The new md_harmonics added, are constructed using only
    /// their period (i.e. all other parametrs are default-constructed).
    /// @param[in] periods A vector<double> containing period values for the
    ///                    harmonics to be added.
    /// @note If you wish to add a harmonic signal in the model, using more info
    /// (e.g. include in- and out-of-phase amplitudes, see the function
    /// add_period).
    void
    add_periods(const std::vector<double>& periods) noexcept
    {
        for (auto i : periods) m_harmonics.emplace_back(i);
    }

    /// This function adds a new md_harmonics instance to this (model) instance.
    /// The new harmonic to be added, can be fully specified.
    /// @param[in] period    The period of the harmonic.
    /// @param[in] in_phase  Amplitude of the in-phase component.
    /// @param[in] out_phase Amplitude of the out-of-phase component.
    /// @param[in] start     Start of the validity interval (default value
    ///                      datetime<T>::min()).
    /// @param[in] stop      End of the validity interval (default value
    ///                      datetime<T>::max()).
    void
    add_period(double period, double in_phase=0e0, double out_of_phase=0e0,
        datetime<T> start=datetime<T>::min(),
        datetime<T> stop=datetime<T>::max())
    noexcept
    {
        m_harmonics.emplace_back(period, start, stop, in_phase, out_of_phase);
    }

    /// Add a jump (i.e. offset) in this model instance. This function will
    /// create a new md_jump instance based on the input parametrs and add it
    /// to this md_model.
    /// @param[in] at  The epoch the jump happened.
    /// @param[in] val The value of the offset (default is 0).
    void
    add_jump(datetime<T> at, double val=0e0) noexcept
    {
        m_jumps.emplace_back(at, val);
    }

    /// Add a velocity change in this model instance. This function will
    /// create a new md_velocity_change instance based on the input parametrs
    /// and add it to this md_model.
    /// @param[in] start   The starting epoch of this velocity change.
    /// @param[in] val     The magnitude of the velocity change.
    /// @param[in] stop    The end epoch of this velocity change; default is
    ///                    ngpt::datetime<T>::max(), i.e. this velocity change
    ///                    starts at epoch start and never stops.
    void
    add_velocity_change(datetime<T> start, double val=0e0,
        ngpt::datetime<T> stop = ngpt::datetime<T>::max())
    noexcept
    {
        m_vel_changes.emplace_back(start, stop, val);
    }

    void
    add_earthquake(datetime<T> start, psd_model m, double a1=1e-3, double t1=1e-1,
        double a2=1e-3, double t2=1e-1)
    noexcept
    {
        m_earthqs.emplace_back(start, m, a1, t1, a2, t2);
    }

    void
    add_earthquake(const md_earthquake<T>& e)
    noexcept
    { m_earthqs.emplace_back(e); }

    /// Remove an earthquake from the earthquake events vector, given a date.
    /// The function will search for an earthquake that has happened at this
    /// date and if found, it will remove it from the model.
    void
    erase_earthquake_at(const datetime<T>& t) noexcept
    {
      for (auto it = m_earthqs.begin(); it != m_earthqs.end(); ++it) {
        if ( it->start() == t ) {
          m_earthqs.erase(it);
          return;
        }
      }
      /* I do not understand why this does not work! WTF???
         auto it = std::find_if(m_earthqs.begin(), m_earthqs.end(),
         [&t](const md_earthquake<T>& a){a.start() == t;});
         if (it!=m_earthqs.end()) m_earthqs.erase(it);
       */
    };
    
    void
    change_earthquake_at(const datetime<T>& t, psd_model md, double a1=0e0, double t1=1.0e0, double a2=0e0, double t2=1.0e0) noexcept
    {
      for (auto it = m_earthqs.begin(); it != m_earthqs.end(); ++it) {
        if ( it->start() == t ) {
          it->set_psd_type(md);
          it->a1() = a1;
          it->a2() = a2;
          it->t1() = t1;
          it->t2() = t2;
          return;
        }
      }
      /* I do not understand why this does not work! WTF???
         auto it = std::find_if(m_earthqs.begin(), m_earthqs.end(),
         [&t](const md_earthquake<T>& a){a.start() == t;});
         if (it!=m_earthqs.end()) m_earthqs.erase(it);
       */
    };

    /// Return the total number of parameters of this instance.
    /// (i.e. 2 for the linear term + 1 for each velocity change + 2 for each
    /// harmonic + 1 for each jump + 4 for each PSD/earthquake)
    std::size_t
    parameters() const noexcept
    {   
        int psd = 0;
        for (auto j = m_earthqs.cbegin(); j != m_earthqs.cend(); ++j)
            psd += j->parameters();

        return  2                        // 2 params for linear terms
            + m_jumps.size()           // 1 param per jump
            + m_vel_changes.size()     // 1 param per vel. change
            + 2*m_harmonics.size()     // 2 params per harmonic
            + psd;
    }

    /// This function is usefull for filling values in A and b matrices, for
    /// least squares parameter estimation (given that we model Ax = b).
    /// Given the A and b matrices and a row number, this function will fill in
    /// the elements of A and b (for the given row), based on the parameters
    /// passed in. Hence, to fully fill in the A and b matrices, you need
    /// recursive call to this function, for all rows.
    /// Let's suppose that the number of observations is N and the number of
    /// parameters is M, then A must be of size (nxM) and b must be of size
    /// (nx1).
    /// This function actually adds a new observation equation to the system
    /// Ax = b.
    ///
    /// @param[in] A     An Eigen::MatrixXd of size (NxM), aka the design
    ///                  matrix. Only row number 'row' is going to be affected
    ///                  in this function.
    /// @param[in] b     An Eigen::VectorXd of size (Nx1), aka the observation
    ///                  matrix. Only row number 'row' is going to be affected
    ///                  in this function.
    /// @param[in] t     The (current) epoch of observation number #row.
    /// @param[in] w     The weight of the observation.
    /// @param[in] dt    The delta-time value, aka t_i - t0 where t0 is the
    ///                  central epoch of computation and t_i is the epoch of
    ///                  the observation number #row (type is double).
    /// @param[in] y     The value of the observation number #row.
    /// @param[in] row   The row (aka index) of the A and b matrices, where we
    ///                  will fill the new observation equation.
    ///
    /// @todo The delta-times are scaled by year. I need to do something
    ///       different with this ... e.g. introduce a scale factor.
    void
    assign_row(Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const datetime<T>& t, double w, double dt, double y, std::size_t row)
    const
    {
        std::size_t col {0};
        //  the observation model value based on approximate values; this is
        //+ used in the non-linear case.
        double l {0e0}; 
        double derivs[4];
        
        // coef for constant (linear) velocity i.e. m/year
        A(row, col) = 1.0e0 * w;
        ++col;
        A(row, col) = w * (dt/365.25e0);
        ++col;
        l += x0() + (dt/365.25e0)*vx();
        
        // Harmonic coefficients for each period ...
        for (auto j = m_harmonics.cbegin(); j != m_harmonics.cend(); ++j) {
            if ( t >= j->start()&& t < j->stop() ) {
                // cosinus or phase
                A(row, col)   = std::cos(j->angular_frequency()*dt) * w;
                l += j->in_phase()*std::cos(j->angular_frequency()*dt);
                // sinus or out-of-phase
                A(row, col+1) = std::sin(j->angular_frequency()*dt) * w;
                l += j->out_of_phase()*std::sin(j->angular_frequency()*dt);
            } else {
                A(row, col)   = 0e0;
                A(row, col+1) = 0e0;
            }
            col+=2;
        }

        // Set up jumps ...
        for (auto j = m_jumps.cbegin(); j != m_jumps.cend(); ++j) {
            if ( t >= j->start() ) {
                A(row, col) = w;
                l += j->value();
            } else {
                A(row, col) = 0e0;
            }
            ++col;
        }
        
        // Set up velocity changes ...
        for (auto j = m_vel_changes.cbegin(); j != m_vel_changes.cend(); ++j) {
            if ( t >= j->start() && t < j->stop() ) {
                auto   dti  = ngpt::delta_date(t, j->start());
                double dtj  = dti.as_mjd();
                A(row, col) = w * (dtj/365.25);
                l += j->value()*(dtj/365.25);
            } else {
                A(row, col) = 0e0;
            }
            ++col;
        }
        
        // Set up earthquake PSD corrections
        for (auto j = m_earthqs.cbegin(); j != m_earthqs.cend(); ++j) {
            auto idrv = j->parameters();
            if ( t >= j->start() ) {
                j->diriv_at(t, derivs[0], derivs[1], derivs[2], derivs[3]);
                for (std::size_t i = 0; i < idrv; i++) A(row, col + i) = w * derivs[i];
                l += j->value_at(t);
            } else {
                for (std::size_t i = 0; i < idrv; i++) A(row, col+i) = 0e0;
            }
            col += idrv;
        }

        // Observation Matrix (vector b)
        b(row) = y;                                   // linear case
        if ( !this->is_linear() ) {
            b(row) = b(row) - l;  // non-linear case
        }
        b(row) *= w;

        // All done
        return;
    }

    /// Get the mean epoch of the model (aka the central epoch of computation).
    /// Const version.
    datetime<T>
    mean_epoch() const noexcept
    { return m_mean_epoch; }
    
    /// Get the mean epoch of the model (akak the central epoch of computation).
    /// Non-Const version.
    datetime<T>&
    mean_epoch() noexcept
    { return m_mean_epoch; }

    /// Get the constant linear term (const version).
    double
    x0() const noexcept { return m_x0; }
    
    /// Get the constant linear term (non-const version).
    double&
    x0() noexcept { return m_x0; }

    /// Get the velocity (const version).
    double
    vx() const noexcept { return m_vx; }
    
    /// Get the velocity (non-const version).
    double&
    vx() noexcept { return m_vx; }

    /// Return a copy of the harmonics vector of the model.
    std::vector<md_harmonics<T>>
    harmonics() const noexcept
    { return m_harmonics; }
    
    /// Return a copy of the harmonics vector of the model.
    std::vector<md_harmonics<T>>&
    harmonics() noexcept
    { return m_harmonics; }

    /// Return a copy of the earthquake events vector of the model.
    std::vector<md_earthquake<T>>
    earthquakes() const noexcept
    { return m_earthqs; }

    /// Clear (i.e remove) all earthquakes
    void
    clear_earthquakes() noexcept
    { m_earthqs.clear(); }

    /// Remove a harmonic term given a period
    int
    remove_harmonic_with_period(double per)
    noexcept
    {
        std::size_t original_size = m_harmonics.size();
        m_harmonics.erase(std::remove(m_harmonics.begin(), m_harmonics.end(),
            [per](const md_harmonics<T>& h){return h.period() == per;}),
            m_harmonics.end());
        return original_size - m_harmonics.size();
    }

    /// Get a const iterator to a harmonic given its period.
    typename std::vector<md_harmonics<T>>::const_iterator
    harmonic_with_period(double per) const noexcept
    {
        return std::find_if(m_harmonics.cbegin(), m_harmonics.cend(),
            [per](const md_harmonics<T>& h){return h.period() == per;});
    }

private:
    double                             m_x0,          ///< const linear term.
                                       m_vx,          ///< linear velocity
                                       m_x0_stddev,
                                       m_vx_stddev;
    std::vector<md_jump<T>>            m_jumps;       ///< vector of jumps
    std::vector<md_harmonics<T>>       m_harmonics;   ///< vector of harmonics
    std::vector<md_velocity_change<T>> m_vel_changes; ///< vector of velocity changes
    std::vector<md_earthquake<T>>      m_earthqs;     ///< vector of earthquakes
    datetime<T>                        m_mean_epoch;  ///< mean epoch (aka central computation epoch)
    
    void
    add_solution_vector(const Eigen::VectorXd& x_estim)
    {
        assert( x_estim.size() >= 2 
                && (int)x_estim.size() == (int)this->parameters() );

        std::size_t idx = 0;
        // std::cout<<"\n\t -- x0 = "<<m_x0<<" + "<<x_estim(idx);
        m_x0 += x_estim(idx); ++idx;
        // std::cout<<"\n\t -- Vx = "<<m_vx<<" + "<<x_estim(idx);
        m_vx += x_estim(idx); ++idx;
        for (auto j = m_harmonics.begin(); j != m_harmonics.end(); ++j) {
            j->in_phase()     += x_estim(idx); ++idx;
            j->out_of_phase() += x_estim(idx); ++idx;
        }
        for (auto j = m_jumps.begin(); j != m_jumps.end(); ++j) {
            j->value() += x_estim(idx); ++idx;
        }
        for (auto j = m_vel_changes.begin(); j != m_vel_changes.end(); ++j) {
            j->value() += x_estim(idx); ++idx;
        }
        for (auto j = m_earthqs.begin(); j != m_earthqs.end(); ++j) {
            std::size_t p = j->parameters();
            //std::cout<<"\n\t -- a1 = "<<j->a1()<<" + "<<x_estim(idx);
            j->a1() += x_estim(idx); ++idx;
            if ( p > 1 ) {
                //std::cout<<"\n\t -- t1 = "<<j->t1()<<" + "<<x_estim(idx);
                j->t1() += x_estim(idx); ++idx;
            }
            if ( p > 3 ) {
                //std::cout<<"\n\t -- a2 = "<<j->a2()<<" + "<<x_estim(idx);
                j->a2() += x_estim(idx); ++idx;
                //std::cout<<"\n\t -- t2 = "<<j->t2()<<" + "<<x_estim(idx);
                j->t2() += x_estim(idx); ++idx;
            }
        }
        return;
    }

}; // end class ts_model

template<class T>
    void
    models_to_json(std::ostream& os, const ts_model<T>& x,
            const ts_model<T>& y, const ts_model<T>& z)
{
    os << "{\n\"model_x\":";
    x.dump_as_json(os);
    os << ",\n\"model_y\":";
    y.dump_as_json(os);
    os << ",\n\"model_z\":";
    z.dump_as_json(os);
    os << "\n}";
}

}// end namespace ngpt

#endif

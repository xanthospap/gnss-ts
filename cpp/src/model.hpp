#ifndef __NGPT_TS_MODEL__
#define __NGPT_TS_MODEL__

// Eigen headers
#include "eigen3/Eigen/Core"

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "event_list.hpp"

namespace ngpt
{

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_jump
{
public:
    explicit
    md_jump(ngpt::datetime<T> start, double offset_in_meters=0) noexcept
    : m_start{start},
      m_offset{offset_in_meters}
    {};

    datetime<T>
    start() const noexcept { return m_start; }

    double
    value() const noexcept { return m_offset; }
    
    double&
    value() noexcept { return m_offset; }

private:
    ngpt::datetime<T> m_start;
    double            m_offset; // meters

}; // md_jump

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_harmonics
{
public:
    explicit
    md_harmonics(double period_in_days, ngpt::datetime<T> start = ngpt::datetime<T>::min(),
            ngpt::datetime<T> stop = ngpt::datetime<T>::max(), double in_phase_val=0,
            double out_phase_val=0)
    noexcept
    : m_start{start},
      m_stop{stop},
      m_afreq{D2PI/period_in_days},
      m_in_phase{in_phase_val},
      m_out_phase{out_phase_val}
    {};


    ngpt::datetime<T>
    start() const noexcept { return m_start; }

    ngpt::datetime<T>
    stop() const noexcept { return m_stop; }

    double
    angular_frequency() const noexcept { return m_afreq; }

    double
    in_phase() const noexcept { return m_in_phase; }

    double
    out_of_phase() const noexcept { return m_out_phase; }
    
    double&
    in_phase() noexcept { return m_in_phase; }

    double&
    out_of_phase() noexcept { return m_out_phase; }

private:
    ngpt::datetime<T> m_start, m_stop;
    double            m_afreq;        // angular frequency i.e. omegas: 2 * pi * frequency)
    double            m_in_phase,     // meters
                      m_out_phase;    // meters

}; // harmonic_coef

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_velocity_change
{
public:
    explicit
    md_velocity_change(ngpt::datetime<T> start, ngpt::datetime<T> stop = ngpt::datetime<T>::max(),
        double new_vel=0)
    noexcept
    : m_start{start},
      m_stop{stop},
      m_newvel{new_vel}
    {};

    ngpt::datetime<T>
    start() const noexcept { return m_start; }

    ngpt::datetime<T>
    stop() const noexcept { return m_stop; }

    double
    value() const noexcept { return m_newvel; }
    
    double&
    value() noexcept { return m_newvel; }
private:

    ngpt::datetime<T> m_start, m_stop;
    double            m_newvel;       // meters/year

}; // velocity_change

template<class T> class ts_model
{
public:

    ts_model() noexcept {};
      
    explicit
    ts_model(const event_list<T>& events) noexcept
    : m_x0{0}, m_vx{0}
    {
        for (auto it = events.it_begin(); it != events.it_end(); ++it)
        {
            if ( it->second == ts_event::jump ) {
                m_jumps.emplace_back(it->first);
            } else if ( it->second == ts_event::velocity_change ) {
                m_vel_changes.emplace_back(it->first);
            } else if ( it->second == ts_event::earthquake ) {
                m_jumps.emplace_back(it->first);
            }
        }
    }

    void
    assign_solution_vector(const Eigen::VectorXd& x_estim)
    {
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
        return;
    }

    void
    add_periods(const std::vector<double>& periods) noexcept
    {
        for (auto i : periods) m_harmonics.push_back(i);
    }

    std::size_t
    parameters() const noexcept
    { return  2 + m_jumps.size() + m_vel_changes.size() + 2*m_harmonics.size(); }

    void
    assign_row(Eigen::MatrixXd& A, Eigen::VectorXd& b, const datetime<T>& current_epoch,
        double weight, double dt, double obs_val, std::size_t row) const
    {
        std::size_t col = 0;
        
        // coef for constant (linear) velocity i.e. m/year
        A(row, col) = 1.0e0 * weight;
        ++col;
        A(row, col) = weight * (dt / 365.25);
        ++col;
        
        // Harmonic coefficients for each period ...
        // typename std::vector<md_harmonics<T>>::const_iterator j;
        for (auto j = m_harmonics.cbegin(); j != m_harmonics.cend(); ++j) {
            // cosinus or phase
            A(row, col) = std::cos(j->angular_frequency() * dt) * weight;
            ++col;
            // sinus or out-of-phase
            A(row, col) = std::sin(j->angular_frequency() * dt) * weight;
            ++col;
        }

        // Set up jumps ...
        for (auto j = m_jumps.cbegin(); j != m_jumps.cend(); ++j) {
            if ( j->start() >= current_epoch ) {
                A(row, col) = weight;
            } else {
                A(row, col) = .0e0;
            }
            ++col;
        }
        
        // Set up velocity changes ...
        for (auto j = m_vel_changes.cbegin(); j != m_vel_changes.cend(); ++j) {
            if ( j->start() >= current_epoch && j->stop() < current_epoch ) {
                A(row, col) = weight * (dt / 365.25);
            } else {
                A(row, col) = .0e0;
            }
            ++col;
        }

        // Observation Matrix (vector b)
        b(row) = obs_val * weight;
        return;
    }

private:
    double                             m_x0,
                                       m_vx;
    std::vector<md_jump<T>>            m_jumps;
    std::vector<md_harmonics<T>>       m_harmonics;
    std::vector<md_velocity_change<T>> m_vel_changes;

}; // end class ts_model


}// end namespace ngpt

#endif

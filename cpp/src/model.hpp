#ifndef __NGPT_TS_MODEL__
#define __NGPT_TS_MODEL__

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "event_list.hpp"

namespace ngpt
{

class ts_model
{
public:

    ts_model() noexcept {};
    
    explicit
    ts_model(const event_list<ts_event>& events) noexcept
    : m_x0{0},
      m_vx{0},
    {
        for (auto it = events.it_begin(); it != events.it_end(); ++it)
        {
            if ( it->second == ts_event::jump ) {
                m_jumps.emplace_back(it->first, 0e0);
            } else if ( it->second == ts_event::velocity_change ) {
                md_velocity_change.emplace_back(it->first, it->first, 0e0);
            } else if ( it->second == ts_event::earthquake ) {
                m_jumps.emplace_back(it->first, 0e0);
            }
        }
    }

    void
    add_periods(const cstd::vector<double>& periods) noexcept
    {
        for (auto i : periods)
            m_harmonics.push_back()
    }

private:
    double                          m_x0, m_vx;
    std::vector<md_jump>            m_jumps;
    std::vector<md_harmonics>       m_harmonics;
    std::vector<md_velocity_change> m_vel_changes;

}; // end class ts_model

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

private:
    ngpt::datetime<T> m_start;
    double            m_offset; // meters

}; // jump

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class md_harmonics
{
public:
    explicit
    md_harmonics(ngpt::datetime<T> start, ngpt::datetime<T> stop,
        double period_in_days, double in_phase_val=0, double out_phase_val=0)
    noexcept
    : m_start{start},
      m_stop{stop},
      m_period{period_in_days},
      m_in_phase{in_phase_val},
      m_out_phase{out_phase_val}
    {};

private:
    ngpt::datetime<T> m_start, m_stop;
    double            m_period;       // e.g. for half a year: 365.25/2
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
    md_velocity_change(ngpt::datetime<T> start, ngpt::datetime<T> stop,
        double new_vel=0)
    noexcept
    : m_start{start},
      m_stop{stop},
      m_newvel{new_vel}
    {};
private:

    ngpt::datetime<T> m_start, m_stop;
    double            m_newvel;       // meters/year

}; // velocity_change

}// end namespace ngpt

#endif

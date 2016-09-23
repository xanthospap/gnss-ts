#ifndef __NGPT_TS_MODEL__
#define __NGPT_TS_MODEL__

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

namespace ngpt
{

class ts_model
{
public:
    ts_model(std::vector<> m_events;)
private:
    double m_x0, m_vx;
    std::vector<jump> m_jumps;
    std::vector<harmonic_coef> m_harmonics;
    std::vector<velocity_changes> m_vel_changes;

}; // end class ts_model

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class jump {
public:
    explicit jump(ngpt::datetime<T> start, double offset_in_meters=0) noexcept
    : m_start{start},
      m_offset{offset_in_meters}
    {};
private:
    ngpt::datetime<T> m_start;
    double m_offset; // meters

}; // jump

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class harmonic_coef {
public:
    explicit harmonic_coef(ngpt::datetime<T> start, ngpt::datetime<T> stop,
        double period_in_days, double in_phase_val=0, double out_phase_val=0) noexcept
    : m_start{start},
      m_stop{stop},
      m_period{period_in_days},
      m_in_phase{in_phase_val},
      m_out_of_phase{out_phase_val}
    {};
private:
    ngpt::datetime<T> m_start, m_stop;
    double m_period; // e.g. for half a year: 365.25/2
    double m_in_phase, m_out_of_phase; // meters

}; // harmonic_coef

template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class velocity_change {
public:
    explicit velocity_change(ngpt::datetime<T> start, ngpt::datetime<T> stop, double new_vel=0) noexcept
    : m_start{start},
      m_stop{stop},
      m_newvel{new_vel}
    {};
private:
    ngpt::datetime<T> m_start, m_stop;
    double m_newvel; // meters/year

}; // velocity_change

}// end namespace ngpt

#endif

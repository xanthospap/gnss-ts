#ifndef __NGPT_TS_MODEL__
#define __NGPT_TS_MODEL__

#include <iostream>

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
    period() const noexcept { return D2PI/m_afreq;}

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
    
    /// \warning earthquake events are translated to jumps. 
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
        ngpt::datetime<T> t0 {modified_julian_day{0}, T{0}};
        m_mean_epoch = t0;
    }

    void
    make_model(datetime<T> start, datetime<T> stop, datetime_interval<T> tinterval,
            std::vector<std::pair<datetime<T>,double>>& results)
    const
    {
        datetime<T> t   = start;
        double value    = 0,
               mean_mjd = m_mean_epoch.as_mjd(),
               dt       = 0;

        while ( t < stop ) {
            dt    = t.as_mjd() - mean_mjd;

            // coef for constant (linear) velocity i.e. m/year
            value = m_x0 + m_vx*(dt/365.25);

            // Harmonic coefficients for each period ...
            for (auto j = m_harmonics.cbegin(); j != m_harmonics.cend(); ++j) {
                if (t >= j->start() && t < j->stop() ) {
                    // cosinus or phase
                    value += j->in_phase()*std::cos(j->angular_frequency()*dt);
                    // sinus or out-of-phase
                    value += j->out_of_phase()*std::sin(j->angular_frequency()*dt);
                }
            }

            // Set up jumps ...
            for (auto j = m_jumps.cbegin(); j != m_jumps.cend(); ++j) {
                if ( t >= j->start() ) {
                    value += j->value();
                }
            }

            // Set up velocity changes ...
            for (auto j = m_vel_changes.cbegin(); j != m_vel_changes.cend(); ++j) {
                if ( t >= j->start() && t < j->stop() ) {
                    value += j->value()*(dt/365.25);
                } 
            }

            t += tinterval;
            results.emplace_back(t, value);
        }
        return;
    }
    
    std::vector<double>
    make_model(const std::vector<datetime<T>>& epochs)
    const
    {
        assert( epochs.size() > 0 );
        double value    = 0,
               mean_mjd = m_mean_epoch.as_mjd(),
               dt       = 0;

        std::vector<double> vals;
        vals.reserve(epochs.size());

        for (const auto& t : epochs) {
            dt    = t.as_mjd() - mean_mjd;

            // coef for constant (linear) velocity i.e. m/year
            value = m_x0 + m_vx*(dt/365.25);

            // Harmonic coefficients for each period ...
            for (auto j = m_harmonics.cbegin(); j != m_harmonics.cend(); ++j) {
                if (t >= j->start() && t < j->stop() ) {
                    // cosinus or phase
                    value += j->in_phase()*std::cos(j->angular_frequency()*dt);
                    // sinus or out-of-phase
                    value += j->out_of_phase()*std::sin(j->angular_frequency()*dt);
                }
            }

            // Set up jumps ...
            for (auto j = m_jumps.cbegin(); j != m_jumps.cend(); ++j) {
                if ( t >= j->start() ) {
                    value += j->value();
                }
            }

            // Set up velocity changes ...
            for (auto j = m_vel_changes.cbegin(); j != m_vel_changes.cend(); ++j) {
                if ( t >= j->start() && t < j->stop() ) {
                    value += j->value()*(dt/365.25);
                } 
            }
            vals.emplace_back(value);
        }
        return vals;
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
    dump(std::ostream& os) const
    {
        if (m_mean_epoch.mjd() != modified_julian_day{0})
            os << "Central Epoch      : " << strftime_ymd_hms<T>(m_mean_epoch)<<"\n";

        os << "Linear Coefficients:";
        os << "\n\t" << m_x0 << "\n\t" << m_vx;
        
        if ( m_harmonics.size() ) os << "\nHarmonic Coefficients:";
        for (auto j = m_harmonics.begin(); j != m_harmonics.end(); ++j) {
            os << "\n\t" << j->in_phase()     << " (in-phase)   ";
            os << "From: " << (j->start()==datetime<T>::min()?"start":strftime_ymd_hms<T>(j->start()) );
            os << " To: "  << (j->stop()==datetime<T>::max()?"end":strftime_ymd_hms<T>(j->stop()) );
            os << " Period: " << j->period();
            os << "\n\t" << j->out_of_phase() << " (out-of-pase)";
        }
        
        if ( m_jumps.size() ) os << "\nJumps/Offsets:";
        for (auto j = m_jumps.begin(); j != m_jumps.end(); ++j) {
            os << "\n\t" << j->value() << " at " << strftime_ymd_hms<T>(j->start());
        }

        if ( m_vel_changes.size() ) os << "\nVelocity Changes:";
        for (auto j = m_vel_changes.begin(); j != m_vel_changes.end(); ++j) {
            os << "\n\t" << j->value();
            os << " From: " << (j->start()==datetime<T>::min()?"start":strftime_ymd_hms<T>(j->start()) );
            os << " To: "  << (j->stop()==datetime<T>::max()?"end":strftime_ymd_hms<T>(j->stop()) );
        }

        return;
    }

    void
    filter_parameters()
    {
        typename std::vector<md_jump<T>>::iterator jit;
        while ( (jit = std::find_if(m_jumps.begin(), m_jumps.end(),
            [](const md_jump<T>& j){ return std::abs(j.value()) < .001;})) != m_jumps.end() )
        {
            m_jumps.erase( jit );
        }
        
        typename std::vector<md_harmonics<T>>::iterator hit;
        while ( (hit = std::find_if(m_harmonics.begin(), m_harmonics.end(),
            [](const md_harmonics<T>& j){ return std::abs(j.in_phase()) < .001 && std::abs(j.out_of_phase()) < .001;})) != m_harmonics.end() )
        {
            m_harmonics.erase( hit );
        }

        return;
    }

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

    void
    add_periods(const std::vector<double>& periods) noexcept
    {
        for (auto i : periods) m_harmonics.emplace_back(i);
    }

    void
    add_period(double period_in_days, datetime<T> start=datetime<T>::min(), datetime<T> stop=datetime<T>::max(),
        double in_phase=0e0, double out_of_phase=0e0) noexcept
    {
        m_harmonics.emplace_back(period_in_days, start, stop, in_phase, out_of_phase);
    }

    void
    add_period(double period_in_days, double in_phase, double out_of_phase) noexcept
    {
        m_harmonics.emplace_back(period_in_days, datetime<T>::min(), datetime<T>::max(), in_phase, out_of_phase);
    }

    void
    add_jump(datetime<T> at, double val=0e0)
    {
        m_jumps.emplace_back(at, val);
    }

    void
    add_velocity_change(datetime<T> start, double val=0e0,
        ngpt::datetime<T> stop = ngpt::datetime<T>::max())
    {
        m_vel_changes.emplace_back(start, stop, val);
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
            if ( current_epoch >= j->start()&& current_epoch < j->stop() ) {
                // cosinus or phase
                A(row, col) = std::cos(j->angular_frequency() * dt) * weight;
                // sinus or out-of-phase
                A(row, col+1) = std::sin(j->angular_frequency() * dt) * weight;
            } else {
                A(row, col)   = 0;
                A(row, col+1) = 0;
            }
            col+=2;
        }

        // Set up jumps ...
        for (auto j = m_jumps.cbegin(); j != m_jumps.cend(); ++j) {
            if ( current_epoch >= j->start() ) {
                A(row, col) = weight;
            } else {
                A(row, col) = .0e0;
            }
            ++col;
        }
        
        // Set up velocity changes ...
        for (auto j = m_vel_changes.cbegin(); j != m_vel_changes.cend(); ++j) {
            if ( current_epoch >= j->start() && current_epoch < j->stop() ) {
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

    datetime<T>
    mean_epoch() const noexcept
    { return m_mean_epoch; }
    
    datetime<T>&
    mean_epoch() noexcept
    { return m_mean_epoch; }

    double
    x0() const noexcept { return m_x0; }
    double&
    x0() noexcept { return m_x0; }
    double
    vx() const noexcept { return m_vx; }
    double&
    vx() noexcept { return m_vx; }

private:
    double                             m_x0,
                                       m_vx;
    std::vector<md_jump<T>>            m_jumps;
    std::vector<md_harmonics<T>>       m_harmonics;
    std::vector<md_velocity_change<T>> m_vel_changes;
    datetime<T>                        m_mean_epoch;

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

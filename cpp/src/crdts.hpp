#ifndef __NGPT_CRD_TIMESERIES__
#define __NGPT_CRD_TIMESERIES__

#include <cmath>
#include "dtcalendar.hpp"
#include "genflags.hpp"
#include "tsflagenum.hpp"
#include "timeseries.hpp"
#include "geodesy.hpp"
#include "car2ell.hpp"

namespace ngpt
{

/// A generic time-series class
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class crdts
{
public:
    /// The specific datetime<T> class we will be using.
    using epoch = ngpt::datetime<T>;
    
    /// Simplify the flag type.
    using tflag = ngpt::flag<ngpt::ts_events>;

    /// Constructor
    explicit crdts(std::string name="") noexcept
    : m_name{name}, m_epochs{}, m_x{}, m_y{}, m_z{}
    {}

    /// Copy constructor.
    crdts(const crdts& ts, std::size_t start=0, std::size_t end=0)
    : m_name{ts.m_name},
      m_epochs{},
      m_x{ts.m_x, start, end},
      m_y{ts.m_y, start, end},
      m_z{ts.m_z, start, end}
    {
        if (!start && !end) {
            m_epochs = ts.m_epochs;
        } else {
            if ( !end ) {
                end = ts.m_epochs.size();
            }
            std::vector<epoch> newvec {ts.m_epochs.cbegin()+start,
                                       ts.m_epochs.cbegin()+end}; 
            m_epochs = std::move(newvec);
        }
        set_epoch_ptr();
    }

    /// Move constructor
    crdts(crdts&& ts)
    : m_name{std::move(ts.m_name)},
      m_epochs{std::move(ts.m_epochs)},
      m_x{std::move(ts.m_x)},
      m_y{std::move(ts.m_y)},
      m_z{std::move(ts.m_z)}
    {
        set_epoch_ptr();
    }

    /// Copy assignment.
    crdts& operator=(const crdts& ts)
    {
        if (this != &ts) {
            m_name = ts.m_name;
            m_epochs = ts.m_epochs;
            m_x = ts.m_x;
            m_y = ts.m_y;
            m_z = ts.m_z;
            set_epoch_ptr();
        }
        return *this;
    }
    
    /// Move assignment.
    crdts& operator=(crdts&& ts)
    {
        if (this != &ts) {
            m_name = std::move(ts.m_name);
            m_epochs = std::move(ts.m_epochs);
            m_x = std::move(ts.m_x);
            m_y = std::move(ts.m_y);
            m_z = std::move(ts.m_z);
            set_epoch_ptr();
        }
        return *this;
    }

    /// Add a crdts data point.
    void add(const epoch& t, double x, double y, double z, double sx=1.0,
        double sy=1.0, double sz=1.0, tflag fx=tflag{}, tflag fy=tflag{},
        tflag fz=tflag{})
    {
        m_epochs.emplace_back(t);
        m_x.add_point(x, sx, fx);
        m_y.add_point(y, sy, fy);
        m_z.add_point(z, sz, fz);
        set_epoch_ptr();
    }

    /// Convert from cartesian to topocentric.
    void cartesian2topocentric() noexcept
    {
        double lat, lon, hgt;
        double cf[9], cf_s[9];
        ngpt::data_point px, py, pz;

        // Reference point is mean value
        ngpt::car2ell(m_x.mean(), m_y.mean(), m_z.mean(), lat, lon, hgt);
        double sinf { std::sin(lat) };
        double cosf { std::cos(lat) };
        double sinl { std::sin(lon) };
        double cosl { std::cos(lon) };
        
        ngpt::detail::car2top_matrix(sinf, sinl, cosf, cosl, cf);
        ngpt::detail::car2top_cov_matrix(sinf*sinf, sinl*sinl, cosf*cosf, cosl*cosl, cf);

#ifdef DEBUG
        assert(
            m_x.size() == m_epochs.size() &&
            m_y.size() == m_epochs.size() &&
            m_z.size() == m_epochs.size() );
#endif
        
        for (std::size_t i = 0; i < m_epochs.size(); ++i) {
            px = m_x[i];
            py = m_y[i];
            pz = m_z[i];
            px.value() = cf[0]*m_x[i].value() + cf[1]*m_y[i].value() + cf[2]*m_z[i].value();
            py.value() = cf[3]*m_x[i].value() + cf[4]*m_y[i].value() + cf[5]*m_z[i].value();
            pz.value() = cf[6]*m_x[i].value() + cf[7]*m_y[i].value() + cf[8]*m_z[i].value();
            double xsigma2 = px.sigma() * px.sigma();
            double ysigma2 = py.sigma() * py.sigma();
            double zsigma2 = pz.sigma() * pz.sigma();
            px.sigma() = std::sqrt(cf[0]*xsigma2 + cf[1]*ysigma2 + cf[2]*zsigma2);
            py.sigma() = std::sqrt(cf[3]*xsigma2 + cf[4]*ysigma2 + cf[5]*zsigma2);
            pz.sigma() = std::sqrt(cf[6]*xsigma2 + cf[7]*ysigma2 + cf[8]*zsigma2);
            m_x[i] = px;
            m_y[i] = py;
            m_z[i] = pz;
        }
        return;
    }

private:

    /// Set the epoch pointer of each timeseries component
    void set_epoch_ptr() noexcept
    {
        m_x.epochs() = &m_epochs;
        m_y.epochs() = &m_epochs;
        m_z.epochs() = &m_epochs;
    }

    std::string        m_name;
    std::vector<epoch> m_epochs;
    timeseries<T>      m_x, m_y, m_z;

}; // end class crdts

} // end namespace ngpt

#endif

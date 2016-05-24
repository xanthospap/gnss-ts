#ifndef __NGPT_CRD_TIMESERIES__
#define __NGPT_CRD_TIMESERIES__

#include "dtcalendar.hpp"
#include "genflags.hpp"
#include "tsflagenum.hpp"
#include "timeseries.hpp"

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
        m_x.add_point(z, sz, fz);
        set_epoch_ptr();
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

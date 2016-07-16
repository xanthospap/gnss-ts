#ifndef __NGPT_CRD_TIMESERIES__
#define __NGPT_CRD_TIMESERIES__

// standard headers
#include <cmath>
#include <algorithm>
#include <initializer_list>

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// ggeodesy headers
#include "ggeodesy/ellipsoid.hpp"
#include "ggeodesy/geodesy.hpp"
#include "ggeodesy/car2ell.hpp"
#include "ggeodesy/vincenty.hpp"
#include "ggeodesy/trnsfdtls.hpp"

// gtms headers
#include "genflags.hpp"
#include "timeseries.hpp"
#include "earthquake_cat.hpp"

namespace ngpt
{

enum class coordinate_type : char
{ cartesian, topocentric, ellipsoidal, unknown };
    
bool is_event(ngpt::flag<ts_event> f) noexcept
{
    return f.check(ts_event::jump)
        || f.check(ts_event::velocity_change)
        || f.check(ts_event::earthquake);
}

std::size_t
split_events(std::vector<ts_event>& events,
    std::initializer_list<ngpt::flag<ts_event>> f) noexcept
{
    events.clear();
    for (auto i = f.begin(); i != f.end(); ++i) {
        if (i->check(ts_event::jump)
            && (std::find(events.cbegin(), events.cend(), ts_event::jump)==events.cend()))
            events.emplace_back(ts_event::jump);
        if (i->check(ts_event::velocity_change)
            && (std::find(events.cbegin(), events.cend(), ts_event::velocity_change)==events.cend()))
            events.emplace_back(ts_event::velocity_change);
        if (i->check(ts_event::earthquake)
            && (std::find(events.cbegin(), events.cend(), ts_event::earthquake)==events.cend()))
            events.emplace_back(ts_event::earthquake);
    }
    return events.size();
}

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
    using tflag = ngpt::flag<ts_event>;

    /// An event is described by the event type and a time-stamp (i.e. epoch).
    using tsevent = std::pair<epoch, ts_event>;

    ///
    using entry = typename timeseries<T, ts_event>::entry;
    
    /// standard ellipsoid
    // using Ell = ngpt::ellipsoid::grs80;

    /// Constructor
    explicit crdts(std::string name="") noexcept
    : m_name{name},
      m_epochs{},
      m_x{}, m_y{}, m_z{},
      m_ctype{coordinate_type::unknown}
    {}

    /// Copy constructor. If needed, the user can specify the start and stop
    /// indexes i.e. construct a crdts from another crdts copying only a portion
    /// of the original crdts.
    /// Only the relevant portion of the events list is going to be copied. The
    /// actual index range copied is [start, end).
    ///
    /// \param[in] start Starting index to start copying from.
    /// \param[in] end   Ending index to stop copying at.
    ///
    crdts(const crdts& ts, std::size_t start=0, std::size_t end=0)
    : m_name{ts.m_name},
      m_epochs{},
      m_x{ts.m_x, start, end},
      m_y{ts.m_y, start, end},
      m_z{ts.m_z, start, end},
      m_events{ts.m_events},
      m_ctype{ts.m_ctype}
    {
        // start and end not given; copy the whole epochs vector
        if (!start && !end) {
            m_epochs = ts.m_epochs;
        // end index not given; end index set to the (original) epochs size.
        } else {
            if ( !end ) { end = ts.m_epochs.size(); }
            std::vector<epoch> newvec {ts.m_epochs.cbegin()+start,
                                       ts.m_epochs.cbegin()+end};
            m_epochs = std::move(newvec);
            // leave the events within the new time interval
            std::vector<tsevent> new_events;
            std::size_t sz = m_epochs.size();
            std::copy_if(m_events.cbegin(), m_events.cend(), new_events.begin(),
                [&](const tsevent& it){it.first>=m_epochs[0] && it.second<=m_epochs[sz-1];});
        }
        set_epoch_ptr();
    }

    /// Move constructor
    crdts(crdts&& ts)
    : m_name{std::move(ts.m_name)},
      m_epochs{std::move(ts.m_epochs)},
      m_x{std::move(ts.m_x)},
      m_y{std::move(ts.m_y)},
      m_z{std::move(ts.m_z)},
      m_events{std::move(ts.m_events)},
      m_ctype{std::move(ts.m_ctype)}
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
            m_events = ts.m_events;
            m_ctype = ts.m_ctype;
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
            m_events = std::move(ts.m_events);
            m_ctype = std::move(ts.m_ctype);
            set_epoch_ptr();
        }
        return *this;
    }


    /// Add a crdts data point.
    void add(const epoch& t, double x, double y, double z, double sx=1.0,
        double sy=1.0, double sz=1.0, tflag fx=tflag{}, tflag fy=tflag{},
        tflag fz=tflag{})
    {
        std::vector<ts_event> events;
        events.reserve(3);
        m_epochs.emplace_back(t);
        m_x.add_point(entry{x, sx, fx});
        m_y.add_point(entry{y, sy, fy});
        m_z.add_point(entry{z, sz, fz});
        if (is_event(fx) || is_event(fy) || is_event(fz)) {
            split_events(events, {fx, fy, fz});
            for (auto i = events.begin(); i != events.end(); ++i) {
                m_events.emplace_back(t, *i);
                // m_events.push_back(tsevent{t,*i});
            }
        }
        set_epoch_ptr();
    }

    /// Return the first date
    datetime<T> first_epoch() const noexcept
    { return m_epochs[0]; }

    /// Return the last date
    datetime<T> last_epoch() const noexcept
    { return m_epochs[m_epochs.size()-1]; }

    /// Return the coordinate type
    coordinate_type crd_type() const noexcept { return m_ctype; }
    
    /// Return the coordinate type
    coordinate_type& crd_type() noexcept { return m_ctype; }

    /// Given an earthquake_catalogue, read it through and apply the earthquakes
    /// of interest. For an earthquake to be applied, the following condition
    /// must be met:
    /// M >= -5.60 + 2.17 * log_10(d), where d is the distance from the station
    /// to the epicenter (in meters). This is taken from \cite{fodits}
    /// The number of applied earthquakes is returned.
    /// To 'apply the earthquakes' means that the instance's events records are
    /// augmented to hold the earthquakes of interest.
    ///
    /// \todo Should i also mark the ts records??
    std::size_t
    apply_earthquake_catalogue(earthquake_catalogue<T>& catalogue)
    {
        catalogue.rewind();
        datetime<T> start {this->first_epoch()},
                    stop {this->last_epoch()};
        earthquake<T> eq;
        std::size_t start_search_at = 0,
                    eq_applied = 0;
        double faz = 0,
               baz = 0;
#ifdef DEBUG
        std::size_t eq_read = 0;
#endif

        double slat, slon, shgt, distance;
        ngpt::car2ell<ngpt::ellipsoid::grs80>(m_x.mean(), m_y.mean(), m_z.mean(),
            slat, slon, shgt);

        while ( catalogue.read_next_earthquake(eq) && eq.epoch() <= stop) {
#ifdef DEBUG
            ++eq_read;
#endif
            if (eq.epoch() >= start) {
                distance = eq.epicenter_distance(slat, slon, faz, baz);
                if ( eq.magnitude() >= -5.6 + 2.17 * std::log10(distance) ) {
                    auto lower = std::lower_bound(m_epochs.begin()+start_search_at,
                                m_epochs.end(), eq.epoch());
                    assert(lower != m_epochs.end());
                    auto index = std::distance(m_epochs.begin(), lower);
                    // m_x.mark(index, ts_events::earthquake);
                    // m_y.mark(index, ts_events::earthquake);
                    // m_z.mark(index, ts_events::earthquake);
                    m_events.emplace_back(eq.epoch(), ts_event::earthquake);
                    ++eq_applied;
#ifdef DEBUG
                    std::cout<<"\tAdding earthquake at "<< eq.epoch().stringify() << " (" <<eq.lat()*180/DPI<<", "<<eq.lon()*180/DPI<<"), of size "<<eq.magnitude()<<"M.\n";
#endif
                }
            }
        }
#ifdef DEBUG
        std::cout<<"\tRead "<<eq_read<<" earthquakes from catalogue\n";
#endif
        sort_event_list();
        return eq_applied;
    }
    
    /// \brief Convert from cartesian to topocentric.
    ///
    /// Usin the mean value as reference point, all points in the time-series
    /// are transformed to topocentric (i.e. vectors from the reference point in
    /// a topocentric reference frame) along with their std. deviations.
    void cartesian2topocentric() noexcept
    {
        double lat, lon, hgt;
        double cf[9], cf_s[9];
        entry px, py, pz;

        // Reference point is mean value
        ngpt::car2ell<ngpt::ellipsoid::grs80>(m_x.mean(), m_y.mean(), m_z.mean(),
            lat, lon, hgt);
        double sinf { std::sin(lat) };
        double cosf { std::cos(lat) };
        double sinl { std::sin(lon) };
        double cosl { std::cos(lon) };
        
        ngpt::detail::car2top_matrix(sinf, sinl, cosf, cosl, cf);
        ngpt::detail::car2top_cov_matrix(sinf*sinf, sinl*sinl, cosf*cosf, cosl*cosl, cf);

#ifdef DEBUG
        assert( m_x.size() == m_y.size() && m_y.size() == m_z.size() );
#endif
        
        for (std::size_t i = 0; i < m_x.size(); ++i) {
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
        m_ctype = coordinate_type::topocentric;
        return;
    }

    /// \brief Sort (in chronological order) the events list and remove duplicates.
    ///
    std::size_t
    sort_event_list() noexcept
    {
        if ( m_events.empty() ) return 0;
        // sort
        std::sort(m_events.begin(), m_events.end(),
            [&](const tsevent& a, const tsevent& b){return a.first < b.first;});
        // remove duplicates
        m_events.erase(std::unique(m_events.begin(), m_events.end()),
            m_events.end());
        return m_events.size();
    }

private:

    /// Set the epoch pointer of each timeseries component.
    void set_epoch_ptr() noexcept
    {
        m_x.epochs() = &m_epochs;
        m_y.epochs() = &m_epochs;
        m_z.epochs() = &m_epochs;
    }

    std::string          m_name;         /// name of the timeseries
    std::vector<epoch>   m_epochs;       /// vector of epochs
    timeseries<T, ngpt::ts_event> m_x, m_y, m_z;  /// the individual components
    std::vector<tsevent> m_events;       /// a vector of events (e.g. earthquakes)
    coordinate_type      m_ctype;        /// the coordinate type

}; // end class crdts

} // end namespace ngpt

#endif

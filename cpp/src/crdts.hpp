#ifndef __NGPT_CRD_TIMESERIES__
#define __NGPT_CRD_TIMESERIES__

// standard headers
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include <fstream>
#include <cstring>
#include <stdexcept>
#ifdef DEBUG
#include <iostream>
#endif

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"

// ggeodesy headers
#include "ggeodesy/ellipsoid.hpp"
#include "ggeodesy/geodesy.hpp"
#include "ggeodesy/car2ell.hpp"
#include "ggeodesy/vincenty.hpp"
#include "ggeodesy/trnsfdtls.hpp"
#ifdef DEBUG
#include "ggeodesy/car2top.hpp"
#endif
#ifdef WMULTIT
#include <thread>
#endif

// gtms headers
#include "genflags.hpp"
#include "timeseries.hpp"
#include "earthquake_cat.hpp"

namespace ngpt
{

enum class coordinate_type : char
{ cartesian, topocentric, ellipsoidal, unknown };
 
// Check if a flag is actually an event   
bool is_event(ngpt::flag<ts_event> f) noexcept
{
    return f.check(ts_event::jump)
        || f.check(ts_event::velocity_change)
        || f.check(ts_event::earthquake);
}

// Split a list of flags to individual events; events are filled in the
// 'events' parameter which is cleared at the start of the function.
// \todo this is kinda stupid; do i really need this ??
std::size_t
split_events(
    std::vector<ts_event>& events,
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
    
    /// Size of the time-series (i.e. number of epochs)
    std::size_t size() const noexcept { return m_epochs.size(); }

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
        sort_events_list();
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
        entry  px, py, pz;
        double dx, dy, dz;

        // Reference point is mean value
        ngpt::car2ell<ngpt::ellipsoid::grs80>(m_x.mean(), m_y.mean(), m_z.mean(),
            lat, lon, hgt);
        double sinf { std::sin(lat) };
        double cosf { std::cos(lat) };
        double sinl { std::sin(lon) };
        double cosl { std::cos(lon) };
        
        ngpt::detail::car2top_matrix(sinf, sinl, cosf, cosl, cf);
        ngpt::detail::car2top_cov_matrix(sinf*sinf, sinl*sinl, cosf*cosf, cosl*cosl, cf_s);

#ifdef DEBUG
        assert( m_x.size() == m_y.size() && m_y.size() == m_z.size() && m_x.size() == size() );
#endif
        
        for (std::size_t i = 0; i < size(); ++i) {
            px = m_x[i];
            py = m_y[i];
            pz = m_z[i];
            dx = m_x.mean() - px.value();
            dy = m_y.mean() - py.value();
            dz = m_z.mean() - pz.value();
            px.value() = -(cf[0]*dx + cf[1]*dy + cf[2]*dz);
            py.value() = -(cf[3]*dx + cf[4]*dy + cf[5]*dz);
            pz.value() = -(cf[6]*dx + cf[7]*dy + cf[8]*dz);
            double xsigma2 = px.sigma() * px.sigma();
            double ysigma2 = py.sigma() * py.sigma();
            double zsigma2 = pz.sigma() * pz.sigma();
            px.sigma() = std::sqrt(cf_s[0]*xsigma2 + cf_s[1]*ysigma2 + cf_s[2]*zsigma2);
            py.sigma() = std::sqrt(cf_s[3]*xsigma2 + cf_s[4]*ysigma2 + cf_s[5]*zsigma2);
            pz.sigma() = std::sqrt(cf_s[6]*xsigma2 + cf_s[7]*ysigma2 + cf_s[8]*zsigma2);
            m_x[i] = px;
            m_y[i] = py;
            m_z[i] = pz;
        }
        m_ctype = coordinate_type::topocentric;
        return;
    }

    bool
    apply_event_list_file(const char* evnf)
    {
        std::ifstream fin (evnf, std::ifstream::in);
        if ( !fin.is_open() ) {
#ifdef DEBUG
            throw std::invalid_argument
            ("Could not open event catalogue file: \""+std::string(evnf)+"\"");
#endif
            return false;
        }

        char line[256];
        char *start;
        epoch t, e_start{this->first_epoch()}, e_stop{this->last_epoch()};
        std::vector<tsevent> events;

        while ( fin.getline(line, 256) )
        {
            if ( *line != '#' && *line != 'Y' && *line != ' ' ) { // if not comment line, empty line, or header ...
                start = line;
                t = ngpt::strptime_ymd_hms<T>(line, &start);
                bool resolved = false;
                for (int i = 0; i < 14 && !resolved; ++i) {
                    if ( *start != ' ' ) {
                        switch (*start) {
                            case 'J' :
                                if (t>=e_start && t<=e_stop) events.emplace_back(t, ts_event::jump);
                                resolved = true;
                                break;
                            case 'j' :
                                if (t>=e_start && t<=e_stop) events.emplace_back(t, ts_event::jump);
                                resolved = true;
                                break;
                            case 'E' :
                                if (t>=e_start && t<=e_stop) events.emplace_back(t, ts_event::earthquake);
                                resolved = true;
                                break;
                            case 'e' :
                                if (t>=e_start && t<=e_stop) events.emplace_back(t, ts_event::earthquake);
                                resolved = true;
                                break;
                            case 'V' :
                                if (t>=e_start && t<=e_stop) events.emplace_back(t, ts_event::velocity_change);
                                resolved = true;
                                break;
                            case 'v' :
                                if (t>=e_start && t<=e_stop) events.emplace_back(t, ts_event::velocity_change);
                                resolved = true;
                                break;
                            default:
#ifdef DEBUG
                                throw std::runtime_error("[ERROR] Invalid event flag in event list file \""+std::string(evnf)+"\""
                                "\nFlag is \""+ (*start) +"\"");
#endif
                                return false;
                        }
                    }
                    ++start;
                }
                if (!resolved) {
#ifdef DEBUG
                    throw std::runtime_error("[ERROR] Invalid line in event list file \""+std::string(evnf)+"\""
                        "\nLine is \""+std::string(line)+"\"");
#endif
                    return false;
                }
                resolved = false;
            }
        }
        if ( !fin.eof() ) {
#ifdef DEBUG
                throw std::runtime_error("Error reading events list file: \""+std::string(evnf)+"\".");
#endif
                return false;
        }

#ifdef DEBUG
    std::cout<<"\nFound and applied "<<events.size()<<" events.";
#endif

        // add the newly created vector to the vector of events
        std::move(events.begin(), events.end(), std::back_inserter(m_events));
        // sort the events list
        sort_events_list();
        // all done!
        return true;
    }

    ///
    auto
    qr_fit(std::vector<double>* i_periods = nullptr)
    {
        // don't let the periods vector be NULL
        std::vector<double> periods;
        if ( i_periods ) periods = *i_periods;

        // sort the events list
        this->sort_events_list();

        // split the events into seperate vectors
        std::vector<epoch> jumps, velchgs, earthqs;
        this->split_events_list(jumps, velchgs, earthqs);

        std::cout<<"\nComponent X:";
        /*auto xv =*/ m_x.qr_ls_solve(jumps, velchgs, periods, 1e-3);
        std::cout<<"\nComponent Y:";
        /* auto yv =*/ m_y.qr_ls_solve(jumps, velchgs, periods, 1e-3);
        std::cout<<"\nComponent Z:";
        /* auto zv =*/ m_z.qr_ls_solve(jumps, velchgs, periods, 1e-3);
    
        /*
        std::vector<double> modelx, modely;
        m_z.make_model_line(first_epoch(), last_epoch(), m_z.central_epoch(),
            zv, jumps, velchgs, *periods, modelx, modely);
        std::ofstream fout ("test.mod");
        for (std::size_t i = 0; i < modelx.size(); ++i)
        fout << modelx[i] << " " << modely[i] << "\n";
        fout.close();
        */

        return;
    }

    // \todo entr is shit just for debuging
    epoch depoch(std::size_t i) const { return m_epochs[i]; }
    std::tuple<entry, entry, entry>
        ddata(std::size_t i) const
    {
        return std::tuple<entry, entry, entry>(m_x[i], m_y[i], m_z[i]);
    }

private:

    /// Set the epoch pointer of each timeseries component.
    void set_epoch_ptr() noexcept
    {
        m_x.epoch_ptr() = &m_epochs;
        m_y.epoch_ptr() = &m_epochs;
        m_z.epoch_ptr() = &m_epochs;
    }
    
    /// \brief Sort (in chronological order) the events list and remove duplicates.
    ///
    std::size_t
    sort_events_list() noexcept
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

    // Split a list of events to individual vectors containing:
    // a: jumps
    // b: velocity changes
    // c: earthquakes
    void
    split_events_list(std::vector<epoch>& jumps, std::vector<epoch>& vel_changes,
    std::vector<epoch>& earthquakes) const noexcept
    {
        jumps.clear();
        vel_changes.clear();
        earthquakes.clear();

        for (auto i = m_events.cbegin(); i != m_events.cend(); ++i) {
            if (i->second == ts_event::jump ) {
                jumps.emplace_back(i->first);
            } else if (i->second == ts_event::velocity_change ) {
                vel_changes.emplace_back(i->first);
            } else if (i->second == ts_event::earthquake ) {
                earthquakes.emplace_back(i->first);
            }
        }
        return;
    }

    std::string          m_name;         /// name of the timeseries
    std::vector<epoch>   m_epochs;       /// vector of epochs
    timeseries<T, ngpt::ts_event> m_x, m_y, m_z;  /// the individual components
    std::vector<tsevent> m_events;       /// a vector of events (e.g. earthquakes)
    coordinate_type      m_ctype;        /// the coordinate type

}; // end class crdts

template<typename T>
    std::ostream& operator<<(std::ostream& os, const crdts<T>& ts)
{
    std::string s;
    std::size_t sz = ts.size();
    for (std::size_t i = 0; i < sz; i++)
    {
        // s = ts.depoch(i).stringify();
        auto t = ts.ddata(i);
        os << ts.depoch(i).as_mjd() << " " << std::get<0>(t).value() << " " << std::get<1>(t).value() << " " << std::get<2>(t).value() << "\n";
    }
    return os;
}

} // end namespace ngpt

#endif

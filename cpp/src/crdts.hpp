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
#include "ggdatetime/datetime_write.hpp"

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
#include "event_list.hpp"
#include "model.hpp"

namespace ngpt
{

enum class coordinate_type : char
{ cartesian, topocentric, ellipsoidal, unknown };
 
/// A generic time-series class
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class crdts
{
public:

    using epoch = typename timeseries<T, pt_marker>::epoch;
    using tflag = typename timeseries<T, pt_marker>::tflag;
    using entry = typename timeseries<T, pt_marker>::entry;

    /// Constructor
    explicit
    crdts(std::string name="") noexcept
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
            std::size_t sz = m_epochs.size();
            m_events = m_events.limit_copy(m_epochs[0], m_epochs[sz-1]);
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
    std::size_t
    size() const noexcept { return m_epochs.size(); }

    /// Add a crdts data point.
    void add(
        const epoch& t,
        double x,          double y,         double z,
        double sx=1.0,     double sy=1.0,    double sz=1.0,
        tflag  fx=tflag{}, tflag fy=tflag{}, tflag fz=tflag{}
        )
    {
        m_epochs.emplace_back(t);
        set_epoch_ptr();
        m_x.add_point(entry{x, sx, fx});
        m_y.add_point(entry{y, sy, fy});
        m_z.add_point(entry{z, sz, fz});
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

    /// Return the event list
    event_list<T>
    events() const noexcept { return m_events; }
    
    /// Return the event list
    event_list<T>&
    events() noexcept { return m_events; }

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
        std::size_t /*start_search_at = 0,*/eq_applied = 0;
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
                    m_events.apply(ts_event::earthquake, eq.epoch());
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
        return eq_applied;
    }
    
    /// \brief Convert from cartesian to topocentric.
    ///
    /// Usin the mean value as reference point, all points in the time-series
    /// are transformed to topocentric (i.e. vectors from the reference point in
    /// a topocentric reference frame) along with their std. deviations.
    void
    cartesian2topocentric() noexcept
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
        assert( m_x.data_pts() == m_y.data_pts() && m_y.data_pts() == m_z.data_pts() && m_x.data_pts() == size() );
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

    void
    apply_event_list_file(const char* evn_file)
    {
        m_events.apply_event_list_file(evn_file, first_epoch(), last_epoch());
    }

    ///
    auto
    qr_fit(ngpt::ts_model<T>& xmodel, ngpt::ts_model<T>& ymodel, ngpt::ts_model<T>& zmodel)
    {
        // a-posteriori std. devs
        double x_stddev, y_stddev, z_stddev;

        std::cout<<"\nComponent X:";
        m_x.qr_ls_solve(xmodel, x_stddev, 1e-3);
        std::cout<<"\nComponent Y:";
        m_y.qr_ls_solve(ymodel, y_stddev, 1e-3);
        std::cout<<"\nComponent Z:";
        m_z.qr_ls_solve(zmodel, y_stddev, 1e-3);
    
        return;
    }

    // \todo entr is shit just for debuging
    epoch depoch(std::size_t i) const { return m_epochs[i]; }
    std::tuple<entry, entry, entry>
        ddata(std::size_t i) const
    {
        return std::tuple<entry, entry, entry>(m_x[i], m_y[i], m_z[i]);
    }
 
#ifdef DEBUG
    // just to show the use of running_window
    void
    test_running_window(datetime_interval<T> window)
    {
        std::size_t sum = 0;
        datetime_interval<T> from, to, vto;

        for (auto rw_it  = m_x.rw_begin(window);
                  !rw_it.hit_the_end();
                  ++rw_it)
        {
            auto first  = rw_it.first(),
                 centre = rw_it.centre(),
                 last   = rw_it.last();

            auto vlast  = rw_it.vlast();

            from = centre.delta_time(first);
            to   = last.delta_time(centre);
            vto  = vlast.delta_time(centre);

            //std::cout << "\nFirst : " << strftime_ymd_hms(first.epoch());
            //std::cout << " Centre: "  << strftime_ymd_hms(centre.epoch());
            //std::cout << " Last : "   << strftime_ymd_hms(vlast.epoch());
            //std::cout<<" [-" << from.days().as_underlying_type() <<", +"<< vto.days().as_underlying_type() <<" (" << to.days().as_underlying_type() << ")]";
            auto data = rw_it.clean_average();
            auto median = rw_it.median();
            auto iqr = rw_it.iqr();
            // std::cout << "\n" << centre.epoch().as_mjd() << " " << data.value() << " " << data.sigma();
            std::cout << " Centre: "  << strftime_ymd_hms(centre.epoch()) << "\n";
            for (auto i = rw_it.first(); i!=rw_it.last(); ++i) {
                std::cout<<i.data().value() << ",";
            }
            std::cout<<"\nMedian: "<<median.value()<<", IQR: "<<iqr.value()<<", Average: "<<data.value()<<"\n";
        }
        return;
    }

    //  just to show the use of timeseries_iterator
    void
    test_iter()
    {
        // auto newts { std::move(m_y) };
        std::size_t sum = 0;
        for (auto it  = m_x.begin();
                  it != m_x.end();
                  ++it)
        {
            ++sum;
            double old_val = it.data().value();
            it.data().value() = -999.9;
            std::cout<<"\nNew element at distance: "<<it.index()<<", old value: "<<old_val<<", new value: "<<it.data().value();
        }
        sum = 0;
        for (auto it  = m_x.begin();
                  it != m_x.end();
                  ++it)
        {
            std::cout<<"\nNew element at distance: "<<it.index()<<", new value: "<<it.data().value();
        }
        std::cout<<"\nFinal sum = "<<sum<<", size = "<<m_x.size();
    }
#endif

    // TODO this should be const!
    std::ostream& dump(std::ostream& os)/*const*/
    {
        auto x_iter = m_x.begin(),
             y_iter = m_y.begin(),
             z_iter = m_z.begin();

        for (; x_iter != m_x.end(); ++x_iter, ++y_iter, ++z_iter) {
            os << x_iter.epoch().as_mjd() << " "
                << x_iter.data().value() << " " 
                << x_iter.data().sigma() << " " 
                << x_iter.data().flag()  << " "
                << y_iter.data().value() << " " 
                << y_iter.data().sigma() << " " 
                << y_iter.data().flag()  << " "
                << z_iter.data().value() << " " 
                << z_iter.data().sigma() << " " 
                << z_iter.data().flag()  << "\n";
        }
        return os;
    }

    std::ostream& dump_model_line(std::ostream& os, const ts_model<T>& modelx,
            const ts_model<T>& modely, const ts_model<T>& modelz)
    {
        datetime_interval<T> dt {ngpt::modified_julian_day{1}, T{0}};
        
        std::vector<std::pair<datetime<T>, double>> v1, v2, v3;
        std::size_t size_hint = last_epoch().as_mjd() - first_epoch().as_mjd();
        v1.reserve(size_hint);
        v2.reserve(size_hint);
        v3.reserve(size_hint);

        modelx.make_model(first_epoch(), last_epoch(), dt, m_x.central_epoch(), v1);
        modely.make_model(first_epoch(), last_epoch(), dt, m_y.central_epoch(), v2);
        modelz.make_model(first_epoch(), last_epoch(), dt, m_z.central_epoch(), v3);

        for (std::size_t i = 0; i < v1.size(); ++i) {
            os << "\n" << v1[i].first.as_mjd() << " " << v1[i].second << " "
                << v2[i].second << " " << v3[i].second;
        }
        return os;
    }

    std::ostream&
    dump_event_list(std::ostream& os) 
    {
        return m_events.dump_event_list(os);
    }

private:

    /// Set the epoch pointer of each timeseries component.
    void set_epoch_ptr() noexcept
    {
        m_x.epoch_ptr() = &m_epochs;
        m_y.epoch_ptr() = &m_epochs;
        m_z.epoch_ptr() = &m_epochs;
    }
    
    std::string          m_name;         /// name of the timeseries
    std::vector<epoch>   m_epochs;       /// vector of epochs
    timeseries<T, ngpt::pt_marker> m_x, m_y, m_z;  /// the individual components
    event_list<T>        m_events;
    coordinate_type      m_ctype;        /// the coordinate type

}; // end class crdts

} // end namespace ngpt

#endif

#ifndef __NGPT_EVENT_LIST__
#define __NGPT_EVENT_LIST__

// c++ standard headers
#include <algorithm>
#include <fstream>
#include <cstring>
#include <stdexcept>

// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"
#include "ggdatetime/datetime_write.hpp"

// gtms headers
#include "genflags.hpp"
#include "earthquake_cat.hpp"

namespace ngpt
{

/// Check if a flag is actually an event, i.e. is any of:
/// - jump,
/// - velocity_change
/// - earthquake
///
bool is_event(ngpt::flag<ts_event> f) noexcept
{
    return f.check(ts_event::jump)
        || f.check(ts_event::velocity_change)
        || f.check(ts_event::earthquake);
}

/// \brief A class to handle time-series events.
///
/// This is a template class, based on (template) parameter:
/// T : which can be any *second type (e.g. second, millisecond, etc) and
///     defines datetime<T>, i.e. datetime resolution, and
///
/// TODO when adding an event (apply()) watch for duplicates
template<class T> class event_list
{
public:

    using epoch   = ngpt::datetime<T>;
    using tflag   = ngpt::flag<ts_event>;
    using event   = std::pair<epoch, ts_event>;
    using earthq  = ngpt::earthquake<T>;

    /// Null constructor.
    event_list() noexcept {};
    
    /// Return a new event_list, which is a calling of the calling instance, but
    /// does not contain events that fall outside the interval [start, stop].
    event_list<T>
    limit_copy(epoch&& start, epoch&& stop) const noexcept
    {
        event_list<T> new_events;
        std::copy_if(m_events.cbegin(), m_events.cend(), new_events.begin(),
                [&](const event& it){it.first>=start && it.second<=stop;});
        return new_events;
    }

    /// Given a flag of type ngpt::flag<ts_event>, apply the individual events
    /// (if any).
    void
    apply(tflag&& f, const epoch& t) noexcept
    {
        if ( is_event(f) ) {
            if ( f.check(ts_event::jump) )
                sorted_insert(ts_event::jump, t);
            if ( f.check(ts_event::velocity_change) )
                sorted_insert(ts_event::velocity_change, t);
            if ( f.check(ts_event::earthquake) )
                sorted_insert(ts_event::earthquake, t);
        }
        return;
    }

    /// Given a ts_event, apply it.
    void
    apply(ts_event f, const epoch& t) noexcept
    {
        switch (f) {
            case ts_event::jump:
                sorted_insert(f, t);
                return;
            case ts_event::velocity_change:
                sorted_insert(f, t);
                return;
            case ts_event::earthquake:
                sorted_insert(f, t);
                return;
        }
        return;
    }
    
    /// Given an event list file, read it through and apply it, i.e. add all
    /// (unique) events to the m_events list.
    void
    apply_event_list_file(const char* evn_file, epoch&& start, epoch&& stop)
    {
        std::ifstream fin (evn_file, std::ifstream::in);
        if ( !fin.is_open() ) {
            throw std::invalid_argument
            ("Could not open event catalogue file: \""+std::string(evn_file)+"\"");
        }

        char  line[256];
        char  *cbegin;
        epoch t;

        while ( fin.getline(line, 256) )
        {
            if ( *line != '#' && *line != 'Y' && *line != ' ' ) { // if not comment line, empty line, or header ...
                cbegin = line;
                t = ngpt::strptime_ymd_hms<T>(line, &cbegin);
                bool resolved = false;
                for (int i = 0; i < 14 && !resolved; ++i) {
                    if ( *cbegin != ' ' ) {
                        switch (*cbegin) {
                            case 'J' :
                                if (t>=start && t<=stop) apply(ts_event::jump, t);
                                resolved = true;
                                break;
                            case 'j' :
                                if (t>=start && t<=stop) apply(ts_event::jump, t);
                                resolved = true;
                                break;
                            case 'E' :
                                if (t>=start && t<=stop) apply(ts_event::earthquake, t);
                                resolved = true;
                                break;
                            case 'e':
                                if (t>=start && t<=stop) apply(ts_event::earthquake, t);
                                resolved = true;
                                break;
                            case 'V' :
                                if (t>=start && t<=stop) apply(ts_event::velocity_change, t);
                                resolved = true;
                                break;
                            case 'v' :
                                if (t>=start && t<=stop) apply(ts_event::velocity_change, t);
                                resolved = true;
                                break;
                            default:
                                throw std::runtime_error("[ERROR] Invalid event flag in event list file \""+std::string(evn_file)+"\""
                                "\nFlag is \""+ (*cbegin) +"\"");
                        }
                    }
                    ++cbegin;
                }
                if ( !resolved ) {
                    throw std::runtime_error("[ERROR] Invalid line in event list file \""+std::string(evn_file)+"\""
                        "\nLine is \""+std::string(line)+"\"");
                }
                resolved = false;
            }
        }
        if ( !fin.eof() ) {
            throw std::runtime_error("Error reading events list file: \""+std::string(evn_file)+"\".");
        }

#ifdef DEBUG
    std::cout<<"\nFound and applied "<<m_events.size()<<" events.";
#endif

        // all done!
        return;
    }

    /// Split the vector of events (aka m_events) to individual vectors per
    /// event (i.e. one vector for jumps, one for velocity changes and one for
    /// outliers.
    /// a: jumps
    /// b: velocity changes
    /// c: earthquakes
    void
    split_event_list(
        std::vector<epoch>& jumps,
        std::vector<epoch>& vel_changes,
        std::vector<epoch>& earthquakes)
    const noexcept
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
    
    /// Write the event list instance to an output stream.
    std::ostream&
    dump_event_list(std::ostream& os) 
    {
        os << "YYYY MM DD HH mm SS **** EVENT *** COMMENT";
        for (auto i = m_events.cbegin(); i != m_events.cend(); ++i) {
            os << "\n" << strftime_ymd_hms(i->first) << "       " 
               << event2char(i->second) << "     ";
        }
        return os;
    }

    typename std::vector<event>::const_iterator
    it_begin() const noexcept
    { return m_events.cbegin(); }
    
    typename std::vector<event>::const_iterator
    it_end() const noexcept
    { return m_events.cend(); }

private:

    /// Insert an event (aka std::pair<datetime<T>, ts_event>) into the
    /// m_events vector in sorted order.
    void
    sorted_insert(ts_event e, const epoch& t)
    {
        event new_event {t, e};       
        m_events.insert(std::upper_bound(m_events.begin(), m_events.end(), new_event,
            [](const event& a, const event& b){return a.first<b.first;}),
            new_event);
        return;
    }

    std::vector<event>   m_events;      /// sorted by epoch
    std::vector<earthq>  m_earthquakes; /// random order

}; // class event_list

} // end namespace ngpt

#endif

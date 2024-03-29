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
#include "tsflag.hpp"
#include "earthquake_cat.hpp"
#include "stalogrd.hpp"

namespace ngpt
{

/// Check if a flag is actually an event, i.e. is any of:
/// - jump,
/// - velocity_change
/// - earthquake
///
bool
is_event(ngpt::flag<ts_event> f) noexcept
{
    return f.check(ts_event::jump)
        || f.check(ts_event::velocity_change)
        || f.check(ts_event::earthquake);
}

template<class T>
    class event
{
public:
    event() noexcept {};

    event(const ngpt::datetime<T>& t, ts_event evn, std::string info=std::string{""})
    noexcept
    : m_epoch(t),
      m_event_type(evn),
      m_info(info)
    {};

    /*
    event(const event<T>& e) noexcept
    : m_epoch(e.epoch()),
      m_event_type(e.event_type()),
      m_info(e.info_str())
    {};

    event(event<T>&& e) noexcept
    : m_epoch(std::move(e.m_epoch)),
      m_event_type(e.m_event_type),
      m_info(std::move(e.m_info))
    {};

    event<T>&
    operator=(const event<T>& e) noexcept
    {
        if (this != &e) {
            m_epoch = e.epoch();
            m_event_type = e.event_type();
            m_info = e.info_str();
        }
        return *this;
    }
    
    event<T>&
    operator=(event<T>&& e) noexcept
    {
        if (this != &e) {
            m_epoch = std::move(e.m_epoch);
            m_event_type = e.event_type();
            m_info = std::move(e.m_info);
        }
        return *this;
    }
    */

    ngpt::datetime<T>
    epoch() const noexcept { return m_epoch; }

    ts_event
    event_type() const noexcept { return m_event_type; }

    std::string
    info_str() const noexcept { return m_info; }

    bool
    operator==(const event<T>& e) const noexcept
    { return (this->m_epoch == e.m_epoch)
          && (this->m_event_type == e.m_event_type); }
    
    bool
    operator!=(const event<T>& e) const noexcept
    { return !(this->operator==(e)); }

private:
    ngpt::datetime<T> m_epoch;
    ts_event m_event_type;
    std::string m_info;
}; // class event

/// A class to handle time-series events (e.g. jumps, earthquakes, etc...).
/// The basic role in this class, is performed by instances of type
/// ngpt::flag<ts_event>, aka tflags. Remember that a tflag has actually on/off
/// values for each of the events/enums in the ts_event enumeration class (so
/// that multiple ts_events can be set at any tflag instance).
///
/// @tparam T Can be any *second type (e.g. second, millisecond, etc) and
///           defines datetime<T>, i.e. datetime resolution, and
///
/// @todo when adding an event (apply()) watch for duplicates
template<class T>
    class event_list
{
public:

    using epoch   = ngpt::datetime<T>;
    using tflag   = ngpt::flag<ts_event>;
    using earthq  = ngpt::earthquake<T>;

    /// Null constructor.
    event_list() noexcept {};
    
    /// Return a new event_list, which is a copy of the calling instance, but
    /// does not contain events that fall outside the interval [start, stop].
    /// @param[in] start The epoch to start copying from
    /// @param[in] stop  The epoch to stop copying at
    /// @return          An event_list<T>, containing all events of the calling
    ///                  instance that have a time-stamp between [start, stop].
    /// @todo            Need to also copy the m_details vector!
    event_list<T>
    limit_copy(const epoch& start, const epoch& stop) const noexcept
    {
        event_list<T> new_events;
        std::copy_if(m_events.cbegin(), m_events.cend(), new_events.begin(),
                [&](const event<T>& it){it.epoch()>=start && it.epoch()<=stop;});
        return new_events;
    }

    /// Given a flag of type ngpt::flag<ts_event>, apply the individual events
    /// (if any). By 'apply', i mean that the individual events (i.e. all enums
    /// set in the flag) are added in chronological order at the instance's
    /// m_events vector.
    /// @param[in] f A tflag (aka ngpt::flag<ts_event>)
    /// @param[in] t An epoch (aka ngpt::datetime<T>), i.e. the date/time when
    ///              the events (i.e. the enums set) happened.
    void
    apply(tflag f, const epoch& t, std::string info=std::string{""}) noexcept
    {
        if ( is_event(f) ) {
            if ( f.check(ts_event::jump) )
                sorted_insert(ts_event::jump, t, info);
            if ( f.check(ts_event::velocity_change) )
                sorted_insert(ts_event::velocity_change, t, info);
            if ( f.check(ts_event::earthquake) )
                sorted_insert(ts_event::earthquake, t, info);
        }
        return;
    }

    /// Given a ts_event, apply it. By 'apply', i mean that the ts_event is
    /// added in chronological order at the instance's m_events vector.
    /// @param[in] f The ts_event to add.
    /// @param[in] t The epoch (aka ngpt::datetime<T>) it happened at.
    void
    apply(const event<T>& e) noexcept
    { this->sorted_insert(e); }

    /// Given a ts_event, apply it. By 'apply', i mean that the ts_event is
    /// added in chronological order at the instance's m_events vector.
    /// @param[in] f The ts_event to add.
    /// @param[in] t The epoch (aka ngpt::datetime<T>) it happened at.
    void
    apply(ts_event f, const epoch& t, std::string info=std::string{""}) noexcept
    {
        this->apply(event<T>{t, f, info});
        /*
        switch (f) {
            case ts_event::jump:
                sorted_insert(f, t, info);
                return;
            case ts_event::velocity_change:
                sorted_insert(f, t, info);
                return;
            case ts_event::earthquake:
                sorted_insert(f, t, info);
                return;
        }
        return;
        */
    }

    /// Given an IGS station log file, read through the receiver and antenna
    /// changes blocks, and apply all subsequent changes.
    ///
    /// @param[in] logfl The name of the logfile.
    void
    apply_stalog_file(const char* igsfl, epoch start, epoch stop)
    {
        ngpt::igs_log log { std::string(igsfl) };
        std::string rec_chng {"Receiver Change"};
        auto rec_changes = log.receiver_changes<T>();
        std::string ant_chng {"Antenna Change"};
        auto ant_changes = log.antenna_changes<T>();
        for (auto& t : rec_changes) {
            if (t >= start && t <= stop) {
                apply(ts_event::jump, t, rec_chng);
            } 
        }
        for (auto& t : ant_changes) {
            if (t >= start && t <= stop) {
                apply(ts_event::jump, t, ant_chng);
            }
        }
        return;
    }
    
    /// Given an event list file, read it through and apply it, i.e. add all
    /// (unique) events to the m_events list.
    /// Structure of the event-list file:
    /// - Each line must have a maximum of 256 characters
    /// - Lines starting with (i.e. the first character is) '#', 'Y' and ' '
    ///   are skipped
    /// - The first field in each (non-skipped) line, must be a date, with the
    ///   format YMD-HMS (see the function ngpt::strptime_ymd_hms<T> for more)
    /// - After the date, and within the next 14 places, the event flag must be
    ///   written, either in upper or in lower case.
    /// Only events within the range [start, stop] are added to the instance.
    /// The events are added in chronological order and duplicates are ignored.
    ///
    /// @param[in] evn_file A (c-) string; the name of the events file (see
    ///                     function description for a valid file format.
    /// @param[in] start    Starting epoch for events we are interested in, i.e.
    ///                     any event recorded before start will not be
    ///                     considered.
    /// @param[in] stop     Ending epoch for events we are interested in, i.e.
    ///                     any event recorded after stop will not be
    ///                     considered.
    void
    apply_event_list_file(const char* evn_file, epoch start, epoch stop)
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
            // if not comment line, empty line, or header ...
            if ( *line != '#' && *line != 'Y' && *line != ' ' ) { 
                cbegin = line;
                t = ngpt::strptime_ymd_hms<T>(line, &cbegin);
                bool resolved = false;
                for (int i = 0; i < 14 && !resolved; ++i) {
                    if ( *cbegin != ' ' ) {
                        switch (*cbegin) {
                            case 'J' :
                                if (t>=start && t<=stop)
                                    apply(ts_event::jump, t);
                                resolved = true;
                                break;
                            case 'j' :
                                if (t>=start && t<=stop)
                                    apply(ts_event::jump, t);
                                resolved = true;
                                break;
                            case 'E' :
                                if (t>=start && t<=stop)
                                    apply(ts_event::earthquake, t);
                                resolved = true;
                                break;
                            case 'e':
                                if (t>=start && t<=stop)
                                    apply(ts_event::earthquake, t);
                                resolved = true;
                                break;
                            case 'V' :
                                if (t>=start && t<=stop)
                                    apply(ts_event::velocity_change, t);
                                resolved = true;
                                break;
                            case 'v' :
                                if (t>=start && t<=stop)
                                    apply(ts_event::velocity_change, t);
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
    /// outliers. Remember than an epoch is actually an ngpt::datetime<T>
    /// instance.
    /// @param[in] jumps       At output, the ordered vector of epochs where
    ///                        jumps happened
    /// @param[in] vel_changes At output, the ordered vector of epochs where
    ///                        velocity_change happened
    /// @param[in] earthquakes At output, the ordered vector of epochs where
    ///                        earthquake happened
    /// @note The input vectors are cleared of all entries, before they are
    ///       filled with the events. If they do hold something at input, it
    ///       will be removed at output.
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
            if (i->event_type() == ts_event::jump ) {
                jumps.emplace_back(i->epoch());
            } else if (i->event_type() == ts_event::velocity_change ) {
                vel_changes.emplace_back(i->epoch());
            } else if (i->event_type() == ts_event::earthquake ) {
                earthquakes.emplace_back(i->epoch());
            }
        }
        return;
    }
    
    /// Write the event list instance to an output stream.
    /// @todo why is this not an '<<' operator ??
    std::ostream&
    dump_event_list(std::ostream& os) 
    {
        os << "YYYY MM DD HH mm SS **** EVENT *** COMMENT";
        for (auto i = m_events.cbegin(); i != m_events.cend(); ++i) {
            os << "\n" << strftime_ymd_hms(i->epoch()) << "       " 
               << event2char(i->event_type()) << "     ";
        }
        return os;
    }

    /// Write the event list instance to an output stream in json format.
    std::ostream&
    dump_event_list_as_json(std::ostream& os)
    const
    {
        std::vector<epoch> jumps, vel_chg, erthq;
        this->split_event_list(jumps, vel_chg, erthq);

        os << "{";
        if ( jumps.size() ) {
            os << "\n\"jumps\": [";
            for (auto i = jumps.cbegin(); i != jumps.cend(); ++i) {
                os << "\n{ \"at\": "<< i->as_mjd() << "}";
                if (i != jumps.cend()-1) os << ",";
            }
            os << "\n]";
        }
        if ( vel_chg.size() ) {
            if (jumps.size()) os <<",";
            os << "\n\"velocity_changes\": [";
            for (auto i = vel_chg.cbegin(); i != vel_chg.cend(); ++i) {
                os << "\n{ \"at\": "<< i->as_mjd() << "}";
                if (i != vel_chg.cend()-1) os << ",";
            }
            os << "\n]";
        }
        if ( erthq.size() ) {
            if ( jumps.size() || vel_chg.size() ) os <<",";
            os << "\n\"earthquakes\": [";
            for (auto i = erthq.cbegin(); i != erthq.cend(); ++i) {
                os << "\n{ \"at\": "<< i->as_mjd() << "}";
                if (i != erthq.cend()-1) os << ",";
            }
            os << "\n]";
        }
        os << "}";

        return os;
    }

    /// Get a const_iterator to the begining of the m_events vector.
    typename std::vector<event<T>>::const_iterator
    it_cbegin() const noexcept
    { return m_events.cbegin(); }
    
    /// Get a const_iterator to end (i.e. one past the end) the of the
    /// m_events vector.
    typename std::vector<event<T>>::const_iterator
    it_cend() const noexcept
    { return m_events.cend(); }
    
    /// Get an iterator to the begining of the m_events vector.
    typename std::vector<event<T>>::const_iterator
    it_begin() noexcept
    { return m_events.begin(); }
    
    /// Get an iterator to end (i.e. one past the end) the of the
    /// m_events vector.
    typename std::vector<event<T>>::iterator
    it_end() noexcept
    { return m_events.end(); }

    /// Size of event list (i.e. number of events)
    std::size_t
    size() const noexcept
    { return m_events.size(); }
    
    /// @brief Replace earthquake sequences with individual earthquakes.
    ///
    /// It may happen (and it actually happens a lot!) that earthquake sequences
    /// occur in high frequency mode. That is, we may have an event list with
    /// a lot of earthquakes in a very limited time span, e.g. 3 to 4 days.
    /// This function will identify these sequencies and replace them with a
    /// single earthquake, the one that has the maximum magnitude.
    ///
    /// @parameter[in] dt A (date)time interval; if two or more earthquakes
    ///                   occur within dt, then they are replaced with a single
    ///                   earthquake, the one with the maximum magnitude.
    void
    filter_earthquake_sequences(const ngpt::datetime_interval<T>& dt)
    {
        if ( !m_events.size() ) return;
        
        // normalize the input datetime_interval (just to be sure).
        auto delta_t = dt;
        delta_t.normalize();

        // find the first earthquake in the list.
        auto it = std::find_if(m_events.begin(), m_events.end(),
            [](const event<T>& e)
            {return e.event_type() == ts_event::earthquake;});
        if ( it == m_events.cend() ) return;

        // iterate through all (valid) earthquakes ...
        while ( it < m_events.end() ) {
            auto t1  = (*it).epoch();
            auto pos = std::distance(m_events.begin(), it);
            // find next earthquake.
            auto e_next = std::find_if(it+1, m_events.end(),
                [](const event<T>& e)
                {return (e.event_type() == ts_event::earthquake);});
            // find next earthquake that is more than dt away.
            auto it_end = std::find_if(it+1, m_events.end(),
                [&delta_t, &t1](const event<T>& e)
                {return (e.event_type() == ts_event::earthquake)
                     && (e.epoch().delta_date(t1) > delta_t);});
            if ( it_end == m_events.end() ) { // have we hit the end?
                break;
            } else if (it_end == e_next) {    // is the next earthquake the next 'valid' one?
                it = e_next;
            } else {
                // we have a sequence between the earthquakes at [it, it_end)
                // find the biggest earthquake within the sequence (max_it).
#ifdef DEBUG
                std::cout<<"\n[DEBUG] Earthquake sequence detected! The following earthquake events:";
#endif
                double max_mag = 0e0;
                auto max_it = it;
                for (auto ij = it; ij != it_end; ++ij) {
                    if ( ij->event_type() == ts_event::earthquake ) {
                        auto eq = resolve_noa_earthquake_catalogue_line<T>(ij->info_str());
#ifdef DEBUG
                        std::cout<<"\n\t\t* "<<ij->info_str();
#endif
                        if (eq.magnitude() > max_mag) {
                            max_mag = eq.magnitude();
                            max_it  = ij;
                        }
                    }
                }
                // store max earthquake (as an event).
#ifdef DEBUG
                std::cout<<"\n\tReplaced with "<<max_it->info_str();
#endif
                event<T> max_event {*max_it};
                // (remove-)erase every earthquake in the sequence.
                std::size_t elements_removed = 0;
                auto pte = std::remove_if(it, it_end,
                    [&elements_removed](event<T> e){
                        if (e.event_type() == ts_event::earthquake) {
                            ++elements_removed;
                            return true;
                        }
                        return false;
                    });
                m_events.erase(pte, pte+elements_removed);
                // (sorted) insert the max earthquake (to replace the sequence).
                this->sorted_insert(max_event);
                // re-establish the iterator.
                it  = m_events.begin()+pos;
#ifdef DEBUG
                auto it__ = std::find_if(m_events.begin()+pos, m_events.end(),
                    [](const event<T>& e)
                    {return (e.event_type() == ts_event::earthquake);});
                assert( it == it__ );
#endif
            } // (else)
        } // (while)
        return;
    }

    /// Just a push_back
    /// Use with extreme care as the event_list must be sorted in chronological
    /// order and this function make absolutely no check!
    void push_back(const event<T>& e) noexcept
    { 
        m_events.push_back(e);
    }

    /// Get the event with index i.
    event<T>
    operator()(std::size_t i) const
    { return m_events[i]; }

    /// Get the event at time t
    event<T>
    event_at(const epoch& t) const
    {
        auto it = std::find_if(m_events.cbegin(), m_events.cend(),
            [t](const event<T>& e){return e.epoch() == t;});
        if (it != m_events.cend()) {
            return *it;
        } else {
            throw std::runtime_error("[ERROR] Not event at that time");
        }
    }

private:

    /// Insert an event into the
    /// m_events vector in sorted (i.e. chronologically) order. If the event
    /// is a duplicate, it will not be added.
    ///
    /// @param[in] e A ts_event to be inserted in the instance (i.e. in the 
    ///              m_events member vector).
    /// @param[in] t The time-stamp of the event (i.e. when it happened).
    /// @return      True if the ts_event is indeed added, or false if it is
    ///              not (because it is a duplicate).
    /// @warning   For this algorith to work properly, the m_events vector must
    ///            be sorted (chronologicaly).
    ///
    bool
    sorted_insert(const event<T>& new_event)
    {
        // easy ... vec is empty, just add the event
        if ( !m_events.size() ) {
            m_events.push_back(new_event);
            return true;
        }

        auto it = std::upper_bound(m_events.begin(), m_events.end(), new_event,
            [](const event<T>& a, const event<T>& b){return a.epoch() < b.epoch();});
        
        // check that this is not a duplicate and insert or push_back
        if ( it == m_events.end() && (*(m_events.end()-1) != new_event) ) {
            m_events.push_back(new_event);
            return true;
        } else {
            if ( *(it-1) != new_event ) {
                 m_events.insert(it, new_event);
                 return true;
            }
        }
        return false;
    }

    /// Insert an event (aka std::pair<datetime<T>, ts_event>) into the
    /// m_events vector in sorted (i.e. chronologically) order. If the event
    /// is a duplicate, it will not be added.
    ///
    /// @param[in] e A ts_event to be inserted in the instance (i.e. in the 
    ///              m_events member vector).
    /// @param[in] t The time-stamp of the event (i.e. when it happened).
    /// @return      True if the ts_event is indeed added, or false if it is
    ///              not (because it is a duplicate).
    /// @warning   For this algorith to work properly, the m_events vector must
    ///            be sorted (chronologicaly).
    ///
    bool
    sorted_insert(ts_event e, const epoch& t, std::string info=std::string(""))
    {
        event<T> new_event {t, e, info};
        return this->sorted_insert(new_event);
    }

    std::vector<event<T>> m_events; ///< Vector of events, sorted by epoch.
}; // class event_list

} // end namespace ngpt

#endif

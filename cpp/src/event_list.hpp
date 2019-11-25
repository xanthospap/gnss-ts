#ifndef __NGPT_EVENT_LIST__
#define __NGPT_EVENT_LIST__

///
/// @file event_list.hpp
///

// c++ standard headers
#include <algorithm>
#include <fstream>
#include <cstring>
#include <stdexcept>
// ggdatetime headers
#include "ggdatetime/dtcalendar.hpp"
#include "ggdatetime/datetime_write.hpp"
// ggeodesy headers
#include "ggeodesy/ellipsoid.hpp"
#include "ggeodesy/geodesy.hpp"
#include "ggeodesy/car2ell.hpp"
// gtms headers
#include "tsflag.hpp"
#include "earthquake_cat.hpp"
#include "stalogrd.hpp"

namespace ngpt
{

/// Check if a flag is actually considered an event.
///
/// @param[in] f  An instance of type ngpt::flag<ts_event>
/// @return True if any of the following ts_event is set:
///     - ts_event::jump,
///     - ts_event::velocity_change
///     - ts_event::earthquake
///     False otherwise.
///
bool
is_event(ngpt::flag<ts_event> f) noexcept
{
  return f.check(ts_event::jump)
    || f.check(ts_event::velocity_change)
    || f.check(ts_event::earthquake);
}

/// @class event
///
/// An event is a .... well an event that has taken place (and possibly affects)
/// a time-series of a certain site. For example, an earthquake that took
/// place in the vicinity of a site, is an event.
/// Such events are examined and taken into account when analyzing times-series.
/// An event is characterized by:
/// * the epoch it happened (as ngpt::datetime<T>)
/// * the event type (as ts_event), and
/// * an optional string containing relevant information
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true.
///
/// @see genflags.hpp
///
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
  class event
{
public:

  /// Default constructor; an empty event
  event() noexcept {};

  /// Constructor
  /// @param[in] t    The epoch the vent happened (type datetime<T>)
  /// @param[in] evn  The type of the event (type ts_event)
  /// @param[in] info Information string (optional). If not gven set to ""
  ///                 (type string)
  event(const ngpt::datetime<T>& t, ts_event evn,
      std::string info=std::string{""})
  noexcept
  : m_epoch(t),
    m_event_type(evn),
    m_info(info)
  {};

  /// Instance of the event
  /// @return The epoch the event took place (type datetime<T>)
  ngpt::datetime<T>
  epoch() const noexcept
  { return m_epoch; }

  /// Event type
  /// @return The event type (type ts_event)
  ts_event
  event_type() const noexcept
  { return m_event_type; }

  /// The information string
  /// @return the information string (type string)
  std::string
  info_str() const noexcept { return m_info; }

  /// Check if two events are the same.
  ///
  /// The function will check if the two events happened at the same instance
  /// and have the same type. **It will not check the information strings**.
  ///
  /// @param[in] e An event instance to check against (type event)
  /// @return True if the two events are of the same type and happened at the
  ///         same instance. False otherwise.
  bool
  operator==(const event<T>& e) const noexcept
  { return (this->m_epoch==e.m_epoch) && (this->m_event_type==e.m_event_type); }
    
  /// Check if two events are different.
  ///
  /// The function will check if the two events happened at the same instance
  /// and have the same type. **It will not check the information strings**.
  ///
  /// @param[in] e An event instance to check against (type event)
  /// @return True if the two events are not of the same type or did not  
  ///         happen at the same instance.
  bool
  operator!=(const event<T>& e) const noexcept
  { return !(this->operator==(e)); }

private:
  ngpt::datetime<T> m_epoch; ///< Instance (epoch) the event took place
  ts_event m_event_type;     ///< Event type
  std::string m_info;        ///< Information string (optional)

}; // class event

/// @class event_list
///
/// A class to handle time-series events (e.g. jumps, earthquakes, etc...).
/// An instance of this type, is actually no more than a vector of event 
/// instances **sorted** by epoch.
///
/// @todo when adding an event (apply()) watch for duplicates
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true.
///
/// @see genflags.hpp
///
template<class T>
  class event_list
{
public:

  /// An epoch is just a short name for ngpt::datetime<T>
  using epoch   = ngpt::datetime<T>;
  using tflag   = ngpt::flag<ts_event>;
  using earthq  = ngpt::earthquake<T>;

  /// Null constructor.
  event_list() noexcept {};

  explicit
  event_list(const std::vector<event<T>>& vec) noexcept
  : m_events(vec)
  {}
    
  /// Return a new event_list, which is a copy of the calling instance, but
  /// does not contain events that fall outside the interval [start, stop].
  /// @param[in] start The epoch to start copying from (inclusive). If
  ///                  datetime<T>::min() is passed in (which is the default
  ///                  value), the copying will start from the begining of
  ///                  the calling instance's events vector.
  /// @param[in] stop  The epoch to stop copying at (inclusive). If
  ///                  datetime<T>::max() is passed in (which is the default
  ///                  value), the copying will stop at the end of the    
  ///                  the calling instance's events vector.
  /// @return          An event_list<T>, containing all events of the calling
  ///                  instance that have a time-stamp between [start, stop].
  event_list<T>
  limit_copy(ngpt::datetime<T>start = ngpt::datetime<T>::min(), 
      ngpt::datetime<T> stop = ngpt::datetime<T>::max())
  const noexcept
  {
    event_list<T> new_events;
    auto sz = m_events.size();
    if (!sz) { return new_events; }

    if (start == ngpt::datetime<T>::min()) start = m_events[0].epoch();
    if (stop  == ngpt::datetime<T>::max()) stop  = m_events[sz-1].epoch();

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
  /// @param[in] info An information string (optional)
  void
  apply(ngpt::flag<ts_event> f, const ngpt::datetime<T>& t, 
      std::string info=std::string{""})
  noexcept
  {
    if (is_event(f)) {
      if (f.check(ts_event::jump))
        sorted_insert(ts_event::jump, t, info);
      if (f.check(ts_event::velocity_change))
        sorted_insert(ts_event::velocity_change, t, info);
      if (f.check(ts_event::earthquake))
        sorted_insert(ts_event::earthquake, t, info);
    }
    return;
  }

  /// Given an event, apply it. By 'apply', i mean that the ts_event is
  /// added in chronological order at the instance's m_events vector.
  /// @param[in] f The event to add.
  void
  apply(const event<T>& e) noexcept
  { this->sorted_insert(e); }

  /// Given a ts_event, apply it. By 'apply', i mean that the ts_event is
  /// added in chronological order at the instance's m_events vector.
  /// @param[in] f The ts_event to add.
  /// @param[in] t The epoch (aka ngpt::datetime<T>) it happened at.
  /// @param[in] info An information string (optional)
  void
  apply(ts_event f, const ngpt::datetime<T>& t,
      std::string info=std::string{""})
  noexcept
  { this->sorted_insert(event<T>{t, f, info}); }

  /// Given an IGS station log file, read through the receiver and antenna
  /// changes blocks, and apply all subsequent changes.
  /// All Receiver and Antenna changes recorded in the log file, will be
  /// added to the instance and marked as "tsevent::jump" events.
  /// Optionally, users can pass in start and stop dates, so only that changes
  /// within this interval are considered.
  ///
  /// @param[in] logfl  The name of the logfile.
  /// @param[in] start  Start of the time interval for considering events
  ///                   (inclusive). Any receiver/antenna change that has
  ///                   happened prior to this date, will **NOT** be
  ///                   considered. Default value is ngpt::datetime<T>::min()
  /// @param[in] stop   Stop of the time interval for considering events
  ///                   (inclusive). Any receiver/antenna change that has
  ///                   happened after this date, will **NOT** be
  ///                   considered. Default value is ngpt::datetime<T>::max()
  void
  apply_stalog_file(const char* igsfl,
      ngpt::datetime<T> start = ngpt::datetime<T>::min(), 
      ngpt::datetime<T> stop  = ngpt::datetime<T>::max(),
      bool apply_rec_changes=false)
  {
    ngpt::igs_log log {std::string(igsfl)};

    if (apply_rec_changes) {
      std::string rec_chng {"Receiver Change"};
      auto rec_changes = log.receiver_changes<T>();
      for (auto& t : rec_changes) {
        if (t>= start && t <= stop) {
          apply(ts_event::jump, t, rec_chng);
        }
      }
    }

    std::string ant_chng {"Antenna Change"};
    auto ant_changes = log.antenna_changes<T>();

    for (auto& t : ant_changes) {
      if (t >= start && t<= stop ) {
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
  ///
  /// @todo add a part of an evn file here !!
  void
  apply_event_list_file(const char* evn_file,
        ngpt::datetime<T> start = ngpt::datetime<T>::min(), 
        ngpt::datetime<T> stop  = ngpt::datetime<T>::max())
  {
    std::ifstream fin (evn_file, std::ifstream::in);
    if (!fin.is_open()) {
      throw std::invalid_argument
        ("[ERROR] apply_event_list_file() :: Could not open event catalogue file: \""+std::string(evn_file)+"\"");
    }

    char  line[256];
    char  *cbegin;
    epoch t;

    while (fin.getline(line, 256)) {
      // if not comment line, empty line, or header ...
      if ( *line != '#' && *line != 'Y' && *line != ' ' ) { 
        cbegin = line;
        t = ngpt::strptime_ymd_hms<T>(line, &cbegin);
        bool resolved = false;
        for (int i=0; i<14 && !resolved; ++i) {
          if (*cbegin != ' ') {
            switch (*cbegin) {
              case 'J' :
              case 'j' :
                if (t>=start && t<=stop) apply(ts_event::jump, t);
                resolved = true;
                break;
              case 'E' :
              case 'e':
                if (t>=start && t<=stop) apply(ts_event::earthquake, t);
                resolved = true;
                break;
              case 'V' :
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
        if (!resolved) {
          throw std::runtime_error("[ERROR] Invalid line in event list file \""+std::string(evn_file)+"\""
              "\nLine is \""+std::string(line)+"\"");
        }
        resolved = false;
      }
    }
    if (!fin.eof()) {
      throw std::runtime_error("Error reading events list file: \""+std::string(evn_file)+"\".");
    }

    // all done!
    return;
  }

  /// Split the vector of events (aka m_events) to individual vectors per
  /// event (i.e. one vector for jumps, one for velocity changes and one for
  /// earthquakes. 
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
        std::vector<ngpt::datetime<T>>& jumps,
        std::vector<ngpt::datetime<T>>& vel_changes,
        std::vector<ngpt::datetime<T>>& earthquakes)
  const noexcept
  {
    jumps.clear();
    vel_changes.clear();
    earthquakes.clear();

    for (auto i=m_events.cbegin(); i!=m_events.cend(); ++i) {
      if (i->event_type() == ts_event::jump) {
        jumps.emplace_back(i->epoch());
      } else if (i->event_type() == ts_event::velocity_change) {
        vel_changes.emplace_back(i->epoch());
      } else if (i->event_type() == ts_event::earthquake) {
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
    for (auto i=m_events.cbegin(); i!=m_events.cend(); ++i) {
      os << "\n" << strftime_ymd_hms(i->epoch()) << "       " 
        << event2char(i->event_type()) << "     "
        << i->info_str();
    }
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
  /// @return Size of individual events in this instance.
  std::size_t
  size() const noexcept
  { return m_events.size(); }

  /// Replace a sequence of earthquakes in the range [start, stop) with the
  /// largest earthquake, that is:
  /// first we identify the largest earthquake (in magnitude) in the range
  /// [start, stop)
  /// then we romove all earthquakes in the range [start, stop) except for the
  /// largest one, and return an iterator to this element (aka the largest
  /// earthquake)
  /// @param[in] Iterator to the first element in range (aka start)
  /// @param[in] Iterator to the one-pass the last element in range (aka stop)
  /// @return    Iterator to the largest earthquake in the (new) event_list; new
  ///            because any earthquake in the range [start, stop) with magnitude
  ///            less that the largest will be deleted.
  /// @note      start must be an iterator to an earthquake event.
  typename std::vector<event<T>>::iterator
  replace_sequence_with_largest(typename std::vector<event<T>>::iterator start,
    typename std::vector<event<T>>::iterator stop)
  {
    auto end = m_events.end();
    if (start==end) return m_events.end();

    // find earthquake with largest magnitude
    typename std::vector<event<T>>::iterator largest_it;
    double max_magnitude = 0e0;
    for (auto it=start; it!=stop; ++it) {
      if (it->event_type() == ts_event::earthquake) {
        auto eq = resolve_noa_earthquake_line<T>(it->info_str().c_str());
#ifdef DEBUG
        //std::cout<<"\n\tRemoving earthquake event: "<<it->info_str();
#endif
        if (eq.magnitude()>max_magnitude) {
          largest_it = it;
          max_magnitude = eq.magnitude();
        }
      }
    }

    // remove all earthquakes but the largest one
    const event<T> max_eq (*largest_it);
    bool max_found = false;
    m_events.erase(std::remove_if(start, stop, [&](const auto& ei){
      if (ei.event_type() == ts_event::earthquake) {
        if (max_eq.info_str() == ei.info_str()) {
          return false;
        } else {
          return true;
        }
      } else {
        return false;
      }}), stop);
    // find the earthquake and get/return it's index/iterator
    largest_it = std::find_if(start, m_events.end(), [&](const auto& ei)
      {return ei.event_type()==ts_event::earthquake && ei.info_str()==max_eq.info_str();}
    );
#ifdef DEBUG
    assert(largest_it->event_type() == ts_event::earthquake);
    assert(largest_it->info_str() == max_eq.info_str());
    //std::cout<<"\n\tInserting earthquake event: "<<largest_it->info_str();
#endif

    return largest_it;
  }

  event_list<T>
  filter_earthquakes(double stax, double stay, double staz, double min_mag=5e0,
    double c1=-5e0, double c2=5.55) const
  {
    if (!m_events.size()) return event_list<T>(*this);

    std::vector<event<T>> vec = this->m_events;
    double slat, slon, shgt,
           distance, faz, baz;
    ngpt::car2ell<ngpt::ellipsoid::grs80>(stax, stay, staz, slat, slon, shgt);

    vec.erase(std::remove_if(vec.begin(), vec.end(), [&](const event<T>& e){
        if (e.event_type()==ts_event::earthquake) {
          auto eq = resolve_noa_earthquake_line<T>(e.info_str().c_str());
          if (eq.magnitude()>=min_mag) {
            distance = eq.epicenter_distance(slat, slon, faz, baz) / 1000e0;
            return !(eq.magnitude() >= c1+c2*std::log10(distance));
          } else {
            return true;
          }
        } else {
          return false;
        }
      }), vec.end());

    return event_list<T>(vec);
  }
    
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
  filter_earthquake_sequences(ngpt::datetime_interval<T> dt)
  {
    if (!m_events.size()) return;

    // normalize the input datetime_interval (just to be sure).
    dt.normalize();

    // find the first earthquake in the list.
    auto it = std::find_if(m_events.begin(), m_events.end(),
        [](const event<T>& e)
        {return e.event_type() == ts_event::earthquake;});
    if (it==m_events.cend()) return;

    // iterate through all (valid) earthquakes ...
    // 'it' is the current earthquake, which happens at epoch 't1', at index
    //+ 'pos'.
    //+ next earthquake is 'e_next'
    while (it < m_events.end()) {
      auto t1  = (*it).epoch();
      auto pos = std::distance(m_events.begin(), it);
      // find next earthquake in list
      auto e_next = std::find_if(it+1, m_events.end(),
          [](const event<T>& e)
          {return (e.event_type() == ts_event::earthquake);}
        );
      if (e_next==m_events.cend()) return;
      // find next earthquake that is more than dt away.
      auto it_end = std::find_if(it+1, m_events.end(),
          [&dt, &t1](const event<T>& e)
          {return (e.event_type() == ts_event::earthquake)
            && (e.epoch().delta_date(t1) > dt);}
        );
      if (it_end == e_next) {    // is the next earthquake the next 'valid' one?
        it = e_next;
      } else {
        // we have a sequence between the earthquakes at [it, it_end)
        // find the biggest earthquake within the sequence (max_it).
#ifdef DEBUG
        //std::cout<<"\n[DEBUG] Earthquake sequence detected! The following earthquake events:";
#endif
        it = replace_sequence_with_largest(it, it_end);
      } // (else)
    } // (while)
    return;
  }

  /*
  /// Get the event with index i.
  event<T>
  operator()(std::size_t i) const
  { return m_events[i]; }*/

  /// Get the event at a given time (if there is any)
  /// @todo this function is stupid and should not exist!
  event<T>
  event_at(const ngpt::datetime<T>& t) const
  {
    auto it = std::find_if(m_events.cbegin(), m_events.cend(),
        [t](const event<T>& e){return e.epoch() == t;});
    if (it != m_events.cend()) {
      return *it;
    } else {
      throw std::runtime_error("[ERROR] Not event at that time");
    }
  }
  
  std::vector<event<T>>&
  get_the_vector() noexcept
  { return this->m_events; }

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
    if (!m_events.size()) {
      m_events.push_back(new_event);
      return true;
    }

    auto it = std::upper_bound(m_events.begin(), m_events.end(), new_event,
        [](const event<T>& a, const event<T>& b){return a.epoch() < b.epoch();});

    // check that this is not a duplicate and insert or push_back
    if (it==m_events.end() && (*(m_events.end()-1)!=new_event) ) {
      m_events.push_back(new_event);
      return true;
    } else {
      if (*(it-1) != new_event) {
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

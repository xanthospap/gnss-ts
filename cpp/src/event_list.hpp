#ifndef __NGPT_EVENT_LIST_HPP__
#define __NGPT_EVENT_LIST_HPP__

///
/// @file event_list.hpp
///

#include "ggdatetime/dtcalendar.hpp"
#include "ts_flag.hpp"
#include <fstream>
#include <vector>

namespace ngpt {

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
class event {
public:
  /// Default constructor; an empty event
  event() noexcept {};

  /// Constructor
  /// @param[in] t    The epoch the vent happened (type datetime<T>)
  /// @param[in] evn  The type of the event (type ts_event)
  /// @param[in] info Information string (optional). If not gven set to ""
  ///                 (type string)
  event(const ngpt::datetime<ngpt::milliseconds> &t, ts_event evn,
        std::string info = std::string{""}) noexcept
      : m_epoch(t), m_event_type(evn), m_info(info){};

  /// Instance of the event
  /// @return The epoch the event took place (type datetime<ngpt::milliseconds>)
  ngpt::datetime<ngpt::milliseconds> epoch() const noexcept { return m_epoch; }

  /// Event type
  /// @return The event type (type ts_event)
  ts_event event_type() const noexcept { return m_event_type; }

  /// The information string
  /// @return the information string (type string)
  std::string info_str() const noexcept { return m_info; }

  /// Check if two events are the same.
  ///
  /// The function will check if the two events happened at the same instance
  /// and have the same type. **It will not check the information strings**.
  ///
  /// @param[in] e An event instance to check against (type event)
  /// @return True if the two events are of the same type and happened at the
  ///         same instance. False otherwise.
  bool operator==(const event &e) const noexcept {
    return (this->m_epoch == e.m_epoch) &&
           (this->m_event_type == e.m_event_type);
  }

  /// Check if two events are different.
  ///
  /// The function will check if the two events happened at the same instance
  /// and have the same type. **It will not check the information strings**.
  ///
  /// @param[in] e An event instance to check against (type event)
  /// @return True if the two events are not of the same type or did not
  ///         happen at the same instance.
  bool operator!=(const event &e) const noexcept {
    return !(this->operator==(e));
  }

  std::string to_string() const noexcept;

private:
  ngpt::datetime<ngpt::milliseconds>
      m_epoch;           ///< Instance (epoch) the event took place
  ts_event m_event_type; ///< Event type
  std::string m_info;    ///< Information string (optional)

}; // event

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
class event_list {
public:
  /// An epoch is just a short name for ngpt::datetime<ngpt::milliseconds>
  // using epoch = ngpt::datetime<ngpt::milliseconds>;
  // using tflag = ngpt::flag<ts_event>;
  // using earthq = ngpt::earthquake<T>;

  /// Null constructor.
  event_list() noexcept {};

  explicit event_list(const std::vector<event> &vec) noexcept : m_events(vec) {}

  /// Return a new event_list, which is a copy of the calling instance, but
  /// does not contain events that fall outside the interval [start, stop].
  /// @param[in] start The epoch to start copying from (inclusive). If
  ///                  datetime<ngpt::milliseconds>::min() is passed in (which
  ///                  is the default value), the copying will start from the
  ///                  begining of the calling instance's events vector.
  /// @param[in] stop  The epoch to stop copying at (inclusive). If
  ///                  datetime<ngpt::milliseconds>::max() is passed in (which
  ///                  is the default value), the copying will stop at the end
  ///                  of the the calling instance's events vector.
  /// @return          An event_list<T>, containing all events of the calling
  ///                  instance that have a time-stamp between [start, stop].
  event_list
  limit_copy(ngpt::datetime<ngpt::milliseconds> start =
                 ngpt::datetime<ngpt::milliseconds>::min(),
             ngpt::datetime<ngpt::milliseconds> stop =
                 ngpt::datetime<ngpt::milliseconds>::max()) const noexcept;

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
  ///                   considered. Default value is
  ///                   ngpt::datetime<ngpt::milliseconds>::min()
  /// @param[in] stop   Stop of the time interval for considering events
  ///                   (inclusive). Any receiver/antenna change that has
  ///                   happened after this date, will **NOT** be
  ///                   considered. Default value is
  ///                   ngpt::datetime<ngpt::milliseconds>::max()
  void apply_stalog_file(const char *igsfl,
                         ngpt::datetime<ngpt::milliseconds> start =
                             ngpt::datetime<ngpt::milliseconds>::min(),
                         ngpt::datetime<ngpt::milliseconds> stop =
                             ngpt::datetime<ngpt::milliseconds>::max(),
                         bool apply_rec_changes = true);

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
  void apply_event_list_file(const char *evn_file,
                             ngpt::datetime<ngpt::milliseconds> start =
                                 ngpt::datetime<ngpt::milliseconds>::min(),
                             ngpt::datetime<ngpt::milliseconds> stop =
                                 ngpt::datetime<ngpt::milliseconds>::max());

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
  split_event_list(std::vector<ngpt::datetime<ngpt::milliseconds>> &jumps,
                   std::vector<ngpt::datetime<ngpt::milliseconds>> &vel_changes,
                   std::vector<ngpt::datetime<ngpt::milliseconds>> &earthquakes)
      const noexcept;

  /// Write the event list instance to an output stream.
  /// @todo why is this not an '<<' operator ??
  std::ostream &dump_event_list(std::ostream &os) const;

  /// Get a const_iterator to the begining of the m_events vector.
  typename std::vector<event>::const_iterator it_cbegin() const noexcept {
    return m_events.cbegin();
  }

  /// Get a const_iterator to end (i.e. one past the end) the of the
  /// m_events vector.
  typename std::vector<event>::const_iterator it_cend() const noexcept {
    return m_events.cend();
  }

  /// Get an iterator to the begining of the m_events vector.
  typename std::vector<event>::const_iterator it_begin() noexcept {
    return m_events.begin();
  }

  /// Get an iterator to end (i.e. one past the end) the of the
  /// m_events vector.
  typename std::vector<event>::iterator it_end() noexcept {
    return m_events.end();
  }

  /// Size of event list (i.e. number of events)
  /// @return Size of individual events in this instance.
  std::size_t size() const noexcept { return m_events.size(); }

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
  ///            because any earthquake in the range [start, stop) with
  ///            magnitude less that the largest will be deleted.
  /// @note      start must be an iterator to an earthquake event.
  /*
  typename std::vector<event>::iterator
  replace_sequence_with_largest(typename std::vector<event>::iterator start,
                                typename std::vector<event>::iterator stop);
  */

  /*
  event_list filter_earthquakes(double stax, double stay, double staz,
                                double min_mag = 5e0, double c1 = -5e0,
                                double c2 = 5.55) const;*/

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
  /*void filter_earthquake_sequences(ngpt::datetime_interval<ngpt::milliseconds>
   * dt);*/

  std::vector<event> &get_the_vector() noexcept { return this->m_events; }

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
  bool sorted_insert(const event &new_event) noexcept;

  /// Insert an event (aka std::pair<datetime<ngpt::milliseconds>, ts_event>)
  /// into the m_events vector in sorted (i.e. chronologically) order. If the
  /// event is a duplicate, it will not be added.
  ///
  /// @param[in] e A ts_event to be inserted in the instance (i.e. in the
  ///              m_events member vector).
  /// @param[in] t The time-stamp of the event (i.e. when it happened).
  /// @return      True if the ts_event is indeed added, or false if it is
  ///              not (because it is a duplicate).
  /// @warning   For this algorith to work properly, the m_events vector must
  ///            be sorted (chronologicaly).
  ///
  bool sorted_insert(ts_event e, const datetime<milliseconds> &t,
                     std::string info = std::string("")) noexcept {
    event new_event{t, e, info};
    return this->sorted_insert(new_event);
  }

  std::vector<event> m_events; ///< Vector of events, sorted by epoch.
};                             // class event_list

} // end namespace ngpt

#endif

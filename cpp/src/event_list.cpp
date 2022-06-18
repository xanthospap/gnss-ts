///
/// @file event_list.hpp
///

#include "event_list.hpp"
#include "earthquake.hpp"
#include "genflags.hpp"
#include "datetime/datetime_read.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "igs_sta_log.hpp"
#include "ts_flag.hpp"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>

// YYYY MM DD HH MM SS EVENT COMMENT
std::string dso::event::to_string() const noexcept {
  std::string str = dso::strftime_ymd_hms(m_epoch, ' ');
  switch (m_event_type) {
  case ts_event::jump:
    str += std::string("   j   ");
    break;
  case ts_event::earthquake:
    str += std::string("   e   ");
    break;
  case ts_event::velocity_change:
    str += std::string("   v   ");
    break;
  default:
    str += std::string("   ?   ");
    break;
  }
  return str + m_info;
}

/// Return a new event_list, which is a copy of the calling instance, but
/// does not contain events that fall outside the interval [start, stop].
dso::event_list dso::event_list::limit_copy(
    dso::datetime<dso::milliseconds> start,
    dso::datetime<dso::milliseconds> stop) const noexcept {
  event_list new_events;
  auto sz = m_events.size();
  if (!sz) {
    return new_events;
  }

  std::copy_if(m_events.cbegin(), m_events.cend(), new_events.m_events.begin(),
               [&start = std::as_const(start),
                &stop = std::as_const(stop)](const event &e) {
                 return e.epoch() >= start && e.epoch() <= stop;
               });
  return new_events;
}

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
///                   dso::datetime<dso::milliseconds>::min()
/// @param[in] stop   Stop of the time interval for considering events
///                   (inclusive). Any receiver/antenna change that has
///                   happened after this date, will **NOT** be
///                   considered. Default value is
///                   dso::datetime<dso::milliseconds>::max()
void dso::event_list::apply_stalog_file(
    const char *igsfl, dso::datetime<dso::milliseconds> start,
    dso::datetime<dso::milliseconds> stop, bool apply_rec_changes) {
  dso::igs_log log{std::string(igsfl)};

  if (apply_rec_changes) {
    std::string rec_chng{"Receiver Change"};
    auto rec_changes = log.receiver_changes();
    for (auto &t : rec_changes) {
      if (t >= start && t <= stop) {
        sorted_insert(ts_event::jump, t, rec_chng);
      }
    }
  }

  std::string ant_chng{"Antenna Change"};
  auto ant_changes = log.antenna_changes();

  for (auto &t : ant_changes) {
    if (t >= start && t <= stop) {
      sorted_insert(ts_event::jump, t, ant_chng);
    }
  }

  return;
}

/// Given an event list file, read it through and apply it, i.e. add all
/// (unique) events to the m_events list.
/// Structure of the event-list file:
/// - Each line must have a maximum of 256 characters
/// - Lines starting with (i.e. the first character is) '#' are skipped
/// - The first field in each (non-skipped) line, must be a date, with the
///   format YMD-HMS (see the function dso::strptime_ymd_hms<T> for more)
/// - After the date, and within the next 7 places, the event flag must be
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
/// This is how an evn file is formated:
/// YYYY MM DD HH MM SS EVENT COMMENT
void dso::event_list::apply_event_list_file(
    const char *evn_file, dso::datetime<dso::milliseconds> start,
    dso::datetime<dso::milliseconds> stop) {
  std::ifstream fin(evn_file, std::ifstream::in);
  if (!fin.is_open()) {
    throw std::invalid_argument("[ERROR] apply_event_list_file() :: Could "
                                "not open event catalogue file: \"" +
                                std::string(evn_file) + "\"");
  }

  char line[256];
  char *cbegin;
  datetime<milliseconds> t;

  if (!fin.getline(line, 256)) {
    std::cerr << "\n[ERROR]@apply_event_list_file() : Failed to read header!";
    throw std::invalid_argument(
        "[ERROR]@apply_event_list_file() : Failed to read header!");
  }

  while (fin.getline(line, 256)) {
    // if not comment line, empty line, or header ...
    if (*line != '#') {
      cbegin = line;
      try {
        t = dso::strptime_ymd_hms<milliseconds>(line, &cbegin);
      } catch (std::invalid_argument &e) {
        std::cerr << "\n[ERROR]@apply_event_list_file() :" << e.what();
        std::cerr << "\n[ERROR] Failed to resolve event list file line: \""
                  << line << "\"";
        std::cerr << "\n[ERROR] Input file: \"" << evn_file << "\"";
        throw e;
      }

      std::string inf_str;
      bool resolved = false;
      for (int i = 0; i < 7 && !resolved; ++i) {
        if (*cbegin != ' ') {
          switch (*cbegin) {
          case 'J':
          case 'j':
            if (t >= start && t <= stop) {
              ++cbegin;
              while (*cbegin && *cbegin == ' ')
                ++cbegin;
              sorted_insert(ts_event::jump, t, std::string(cbegin));
            }
            resolved = true;
            break;
          case 'E':
          case 'e':
            if (t >= start && t <= stop) {
              while (*cbegin && *cbegin == ' ')
                ++cbegin;
              sorted_insert(ts_event::earthquake, t, std::string(cbegin));
            }
            resolved = true;
            break;
          case 'V':
          case 'v':
            if (t >= start && t <= stop) {
              while (*cbegin && *cbegin == ' ')
                ++cbegin;
              sorted_insert(ts_event::velocity_change, t, std::string(cbegin));
            }
            resolved = true;
            break;
          default:
            throw std::runtime_error(
                "[ERROR] Invalid event flag in event list file \"" +
                std::string(evn_file) +
                "\""
                "\nFlag is \"" +
                (*cbegin) + "\"");
          }
        }
        ++cbegin;
      }
      if (!resolved) {
        throw std::runtime_error("[ERROR] Invalid line in event list file \"" +
                                 std::string(evn_file) +
                                 "\""
                                 "\nLine is \"" +
                                 std::string(line) + "\"");
      }
      resolved = false;
    }
  }
  if (!fin.eof()) {
    throw std::runtime_error("Error reading events list file: \"" +
                             std::string(evn_file) + "\".");
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
void dso::event_list::split_event_list(
    std::vector<dso::datetime<dso::milliseconds>> &jumps,
    std::vector<dso::datetime<dso::milliseconds>> &vel_changes,
    std::vector<dso::datetime<dso::milliseconds>> &earthquakes)
    const noexcept {
  jumps.clear();
  vel_changes.clear();
  earthquakes.clear();

  for (auto i = m_events.cbegin(); i != m_events.cend(); ++i) {
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
std::ostream &dso::event_list::dump_event_list(std::ostream &os) const {
  os << "YYYY MM DD HH MM SS EVENT COMMENT\n";
  for (auto it = m_events.cbegin(); it != m_events.cend(); ++it) {
    os << it->to_string() << "\n";
  }
  return os;
}

/*
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
typename std::vector<event>::iterator
replace_sequence_with_largest(typename std::vector<event>::iterator start,
                              typename std::vector<event>::iterator stop) {
  std::cout << "\n[DEBUG] calling function replace_sequence_with_largest()";
  auto end = m_events.end();
  if (start == end)
    return m_events.end();

  // find earthquake with largest magnitude
  typename std::vector<event>::iterator largest_it;
  double max_magnitude = 0e0;
  for (auto it = start; it != stop; ++it) {
    if (it->event_type() == ts_event::earthquake) {
      auto eq = resolve_noa_earthquake_line(it->info_str().c_str());
#ifdef DEBUG
      // std::cout<<"\n\tRemoving earthquake event: "<<it->info_str();
#endif
      if (eq.magnitude() > max_magnitude) {
        largest_it = it;
        max_magnitude = eq.magnitude();
      }
    }
  }

  // remove all earthquakes but the largest one
  const event max_eq(*largest_it);
  bool max_found = false;
  m_events.erase(std::remove_if(start, stop,
                                [&](const auto &ei) {
                                  if (ei.event_type() == ts_event::earthquake) {
                                    if (max_eq.info_str() == ei.info_str()) {
                                      return false;
                                    } else {
                                      return true;
                                    }
                                  } else {
                                    return false;
                                  }
                                }),
                 stop);
  // find the earthquake and get/return it's index/iterator
  largest_it = std::find_if(start, m_events.end(), [&](const auto &ei) {
    return ei.event_type() == ts_event::earthquake &&
           ei.info_str() == max_eq.info_str();
  });
#ifdef DEBUG
  assert(largest_it->event_type() == ts_event::earthquake);
  assert(largest_it->info_str() == max_eq.info_str());
  // std::cout<<"\n\tInserting earthquake event: "<<largest_it->info_str();
#endif

  return largest_it;
}
*/

/*
event_list filter_earthquakes(double stax, double stay, double staz,
                              double min_mag = 5e0, double c1 = -5e0,
                              double c2 = 5.55) const {
  std::cout << "\n[DEBUG] calling function filter_earthquakes()";
  if (!m_events.size())
    return event_list(*this);

  std::vector<event> vec = this->m_events;
  double slat, slon, shgt, distance, faz, baz;
  dso::car2ell<dso::ellipsoid::grs80>(stax, stay, staz, slat, slon, shgt);

  vec.erase(
      std::remove_if(
          vec.begin(), vec.end(),
          [&](const event &e) {
            if (e.event_type() == ts_event::earthquake) {
              auto eq = resolve_noa_earthquake_line(e.info_str().c_str());
              if (eq.magnitude() >= min_mag) {
                distance = eq.epicenter_distance(slat, slon, faz, baz) / 1000e0;
                return !(eq.magnitude() >= c1 + c2 * std::log10(distance));
              } else {
                return true;
              }
            } else {
              return false;
            }
          }),
      vec.end());

  return event_list(vec);
}
*/

// TODO this needs to be re-implemented
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
/*
void dso::event_list::filter_earthquake_sequences(
    dso::datetime_interval<dso::milliseconds> dt) {
  if (!m_events.size())
    return;

  // normalize the input datetime_interval (just to be sure).
  dt.normalize();

  // find the first earthquake in the list.
  auto it = std::find_if(m_events.begin(), m_events.end(), [](const event &e) {
    return e.event_type() == ts_event::earthquake;
  });
  if (it == m_events.cend())
    return;

  // iterate through all (valid) earthquakes ...
  // 'it' is the current earthquake, which happens at epoch 't1', at index
  //+ 'pos'.
  //+ next earthquake is 'e_next'
  while (it < m_events.end()) {
    auto t1 = (*it).epoch();
    auto pos = std::distance(m_events.begin(), it);
    // find next earthquake in list
    auto e_next = std::find_if(it + 1, m_events.end(), [](const event &e) {
      return (e.event_type() == ts_event::earthquake);
    });
    if (e_next == m_events.cend())
      return;
    // find next earthquake that is more than dt away.
    auto it_end =
        std::find_if(it + 1, m_events.end(), [&dt, &t1](const event &e) {
          return (e.event_type() == ts_event::earthquake) &&
                 (e.epoch().delta_date(t1) > dt);
        });
    if (it_end == e_next) { // is the next earthquake the next 'valid' one?
      it = e_next;
    } else {
      // we have a sequence between the earthquakes at [it, it_end)
      // find the biggest earthquake within the sequence (max_it).
#ifdef DEBUG
      // std::cout<<"\n[DEBUG] Earthquake sequence detected! The following
      // earthquake events:";
#endif
      it = replace_sequence_with_largest(it, it_end);
    } // (else)
  }   // (while)
  return;
}
*/

/// Insert an event into the m_events vector in sorted (i.e. chronologically)
// order. If the event is a duplicate, it will not be added.
bool dso::event_list::sorted_insert(const event &new_event) noexcept {
  // easy ... vec is empty, just add the event
  if (!m_events.size()) {
    m_events.emplace_back(new_event);
    return true;
  }

  auto it = std::lower_bound(
      m_events.begin(), m_events.end(), new_event,
      [](const event &a, const event &b) { return a.epoch() < b.epoch(); });

  if (it == m_events.end()) {
    m_events.emplace_back(new_event);
    return true;
  } else {
    if (*it != new_event) {
      m_events.insert(it, new_event);
      return true;
    } else {
#ifdef TS_DEBUG
      std::cerr << "\n[DEBUG] Event not added to list cause its a duplictate!";
#endif
    }
  }

  return false;
}

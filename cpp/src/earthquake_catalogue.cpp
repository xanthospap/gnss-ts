#include "earthquake_catalogue.hpp"
#include "datetime/datetime_read.hpp"
#include "geodesy/units.hpp"
#include <exception>
#include <iostream>

/// @brief Constructor using a filename.
dso::earthquake_catalogue::earthquake_catalogue(const std::string &filename)
    : m_filename{filename}, m_ifs{filename.c_str(), std::ifstream::in},
      m_end_of_header{0} {
  if (!m_ifs.is_open()) {
    throw std::invalid_argument(
        "earthquake_catalogue: Could not open earthquake catalogue file: \"" +
        filename + "\"");
  }
  if (!read_header()) {
    m_ifs.close();
    throw std::invalid_argument("earthquake_catalogue: Could not read header "
                                "in earthquake catalogue file: \"" +
                                filename + "\"");
  }
}

/// @brief Resolve a NOA earthquake catalogue file line.
///
/// This function will resolve a NOA earthquake catalogue file line and
/// return the individual fields. The line follows the format:
/// YYYY OOO DD   HH MM SS.S   LAT     LON     DEPTH(km)  MAGNITUDE
/// where 'OOO' is the month as 3-char uppercase string, LAT and LON are
/// given in decimal degrees with a precision of e-2, depth is given in
/// integer km and magnitude in M with precision 1e-1. Example:
/// 1964 FEB 24   23 30 25.0   38.90   23.90   10         5.3
///
/// @parameter[in] line A line of the earthquake catalogue file (to be
///                     resolved (char*).
/// @return             An earthquake instance, comprised of the information
///                     resolved (aka earthquake<T>).
///
/// @note The information are transformed to the 'correct units' for the
///       instance. That is, degrees are transformed to radians, km to
///       meters.
dso::earthquake resolve_noa_earthquake_line(const char *line) {
  float info[4];
  const char *start(line);
  char *end;
  dso::datetime<dso::milliseconds> eph;
  try {
    eph = dso::strptime_yod_hms<dso::milliseconds>(line, &end);
  } catch (std::invalid_argument &e) {
    std::cerr << "\n[ERROR]@resolve_noa_earthquake_line() : " << e.what();
    std::cerr << "\n[ERROR] Invalid date format at line: \"" << line << "\"\n";
    throw e;
  }
  start = end;

  for (int i = 0; i < 4; ++i) {
    info[i] = std::strtod(start, &end);
    if (errno == ERANGE || start == end) {
      errno = 0;
      throw std::invalid_argument(
          "resolve_noa_earthquake_line: Invalid line: \"" + std::string(line) +
          "\" (argument #" + std::to_string(i + 1) + ") in catalogue file.");
    }
    start = end;
  }

  dso::earthquake eqt{eph, dso::deg2rad<double>(info[0]),
                       dso::deg2rad<double>(info[1]), info[2] / 1000e0,
                       info[3]};
  return eqt;
}

/// @brief Read and return the next earthquake
int dso::earthquake_catalogue::read_next_earthquake(earthquake &eq) noexcept {
  char line[earthquake_catalogue_detail::MAX_CHARS_IN_LINE];

  if (!m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE)) {
    if (!m_ifs.eof()) {
      return -1;
    }
    return 1;
  }

  try {
    eq = resolve_noa_earthquake_line(line);
  } catch (std::exception &e) {
    return -2;
  }
  return 0;
}

/// Navigate to a position in the catalogue file, such that the next
/// earthquake to be read in will be at an epoch >= to the given epoch
int dso::earthquake_catalogue::goto_epoch(
    dso::datetime<dso::milliseconds> eph) noexcept {
  dso::earthquake eq;
  std::ifstream::pos_type _pos = m_end_of_header;
  if (m_ifs.tellg() != m_end_of_header)
    this->rewind();
  int status = 0;
  while (!(status = read_next_earthquake(eq))) {
    if (eq.m_epoch >= eph)
      break;
    else
      _pos = m_ifs.tellg();
  }
  if (status != 0)
    return status;
  m_ifs.seekg(_pos, std::ios::beg);
  return 0;
}

/// @brief Read header records.
bool dso::earthquake_catalogue::read_header() noexcept {
  char line[earthquake_catalogue_detail::MAX_CHARS_IN_LINE];
  const char *line1 =
      "      DATE         TIME     LAT.   LONG.  DEPTH    MAGNITUDE";
  const char *line2 =
      "                   (GMT)    (N)    (E)    (km)       (Local)";

  if (!m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE) ||
      std::strncmp(line, line1, std::strlen(line1))) {
    return false;
  }
  if (!m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE) ||
      std::strncmp(line, line2, std::strlen(line2))) {
    return false;
  }
  m_end_of_header = m_ifs.tellg();
  return true;
}

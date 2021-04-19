#include "igs_sta_log.hpp"
#include "ggdatetime/datetime_read.hpp"
#include "ggdatetime/datetime_write.hpp"
#include <cassert>
#include <ggdatetime/dtfund.hpp>

/// Check if a string is empty (aka full of whitespace chars).
///
/// @param[in] str A C-string
/// @return    true if string is empty, false otherwise.
int str_is_empty(const char *str) noexcept {
  while (*str) {
    if (!isspace((unsigned char)*str))
      return false;
    ++str;
  }
  return true;
}

/// Resolve a datetime from a string of type: 'CCYY-MM-DDThh:mmZ' or
/// 'CCYY-MM-DD'
///
/// @param[in] str A c-string containing a date of type 'CCYY-MM-DDThh:mmZ' or
///                'CCYY-MM-DD'.
/// @return The datetime represented by the input string.
ngpt::datetime<ngpt::milliseconds> strptime_log(const char *str,
                                                char **stop = nullptr) {
  errno = 0;
  char *end;
  const char *start = str;
  int ints[5], sz = strlen(str);

  for (int i = 0; i < 5; ++i) {
    ints[i] = static_cast<int>(std::abs(std::strtol(start, &end, 10)));
    if (end - str >= sz && i == 3) { /* CCYY-MM-DD */
      ints[3] = 0;
      ints[4] = 0;
      i = 10;
    }
    if (errno == ERANGE) {
      errno = 0;
      throw std::invalid_argument(
          "igs_log_details::strptime_log: Invalid date format: \"" +
          std::string(str) + "\" (argument #" + std::to_string(i + 1) + ").");
    }
    start = end + 1;
  }

  if (stop)
    *stop = end - 1;
  return ngpt::datetime<ngpt::milliseconds>{
      ngpt::year{ints[0]},  ngpt::month{ints[1]},   ngpt::day_of_month{ints[2]},
      ngpt::hours{ints[3]}, ngpt::minutes{ints[4]}, ngpt::milliseconds{0}};
}

ngpt::igs_log::igs_log(const std::string &filename)
    : m_filename{filename}, m_ifs{filename.c_str(), std::ifstream::in} {
  if (!m_ifs.is_open()) {
    throw std::invalid_argument("igs_log: Could not open log file: \"" +
                                filename + "\"");
  }
}

std::vector<ngpt::datetime<ngpt::milliseconds>>
ngpt::igs_log::receiver_changes() {
  std::size_t dummy{0};
  constexpr std::size_t MAX_CHARS{256}, MAX_LINES{10000};
  const char *receiver_info_block = "3.   GNSS Receiver Information";

  char line[MAX_CHARS];
  ngpt::datetime<ngpt::milliseconds> from,
      to /*,
prev_from = ngpt::datetime<ngpt::milliseconds>::max()*/
      ;
  std::string rectype;
  std::vector<datetime<ngpt::milliseconds>> vdt;

  // go to the begining of the file.
  this->rewind();

  // read all lines untill we reach the block: '3.   GNSS Receiver
  // Information'
  while (m_ifs.getline(line, MAX_CHARS) && strcmp(receiver_info_block, line)) {
    if (++dummy >= MAX_LINES) {
      throw std::invalid_argument("receiver_changes: Cannot find block "
                                  "\'GNSS Receiver Information\'!");
    }
  }
  // confirm that we are on the right line.
  if (strcmp(receiver_info_block, line)) {
    throw std::invalid_argument(
        "receiver_changes: Failed to find receiver info block in log file.");
  }
  // keep on reading receiver blocks ...
  int answer;
  while (!(answer = this->read_receiver_block(rectype, from, to))) {
    assert(to > from /*&& to <= prev_from*/);
    // prev_from = from;
    if (to != ngpt::datetime<ngpt::milliseconds>::max()) {
      vdt.push_back(to);
    }
  }
  // hopefully no error encountered while resolving blocks
  if (answer < 0) {
    std::cerr << "\n[ERROR] Got back " << answer;
  }

  // all done
  return vdt;
}

std::vector<ngpt::datetime<ngpt::milliseconds>>
ngpt::igs_log::antenna_changes() {
  std::size_t dummy{0};
  constexpr std::size_t MAX_CHARS{256}, MAX_LINES{10000};
  const char *antenna_info_block = "4.   GNSS Antenna Information";

  char line[MAX_CHARS];
  ngpt::datetime<ngpt::milliseconds> from, to;
  std::string anttype;
  std::vector<datetime<ngpt::milliseconds>> vdt;

  // go to the begining of the file.
  this->rewind();

  // read all lines untill we reach the block: '4.   GNSS Antenna Information'
  while (m_ifs.getline(line, MAX_CHARS) && strcmp(antenna_info_block, line)) {
    if (++dummy >= MAX_LINES) {
      throw std::invalid_argument(
          "antenna_changes: Cannot find block \'GNSS Antenna Information\'!");
    }
  }
  // confirm that we are on the right line.
  if (strcmp(antenna_info_block, line)) {
    throw std::invalid_argument(
        "antenna_changes: Failed to find antenna info block in log file.");
  }
  // keep on reading antenna blocks ...
  int answer;
  while (!(answer = this->read_antenna_block(anttype, from, to))) {
    assert(to > from /*&& to <= prev_from*/);
    if (to != ngpt::datetime<ngpt::milliseconds>::max()) {
      vdt.push_back(to);
    }
  }
  // hopefully no error encountered while resolving blocks
  if (answer < 0) {
    std::cerr << "\n[ERROR] Got back " << answer;
  }

  // all done
  return vdt;
}

/// Resolve an igs station log's file receiver block (aka the
/// '3.   GNSS Receiver Information' block). This function expects that the
/// instance's stream is ready to read a line of type:
/// '3.x  Receiver Type            : (A20, from rcvr_ant.tab; see
/// instructions)' If the first line read is something different, an error is
/// signaled. It then will read all recorded information and return the
/// important fields. The block should be formatted exactly as described at
/// the IGS specifications (see
/// ftp://igs.org/pub/station/general/sitelog_instr.txt). The function (if
/// successeful), will leave the stream after having read all info of the
/// block, including any 'Additional Information'. That means that the next
/// line to be read, should be an empty line.
///
/// @param[out] receiver_type The receiver type specified in this block (line
///                           'Receiver Type')
/// @param[out] start         Start time of this block's validity (line
///                           'Date Installed'). If this is set to the string
///                           'CCYY-MM-DDThh:mmZ' inside the block, then
///                           start is set to the min
///                           datetime<ngpt::milliseconds>. The same goes for an
///                           empty string.
/// @param[out] stop          Stop time of this block's validity (line
///                           'Date Removed). If this is set to the string
///                           'CCYY-MM-DDThh:mmZ' inside the block, then
///                           stop is set to the max
///                           datetime<ngpt::milliseconds>. The same goes for an
///                           empty string.
/// @return  The function returns an integer, signaling the following:
///           - int < 0 error
///           - int = 0 all ok
///           - int > 0 block 3.x (skipped; no more blocks to read)
int ngpt::igs_log::read_receiver_block(
    std::string &receiver_type, ngpt::datetime<ngpt::milliseconds> &start,
    ngpt::datetime<ngpt::milliseconds> &stop) {
  constexpr std::size_t MAX_CHARS{256}, MAX_LINES{1000};
  char line[MAX_CHARS];

  m_ifs.getline(line, MAX_CHARS); // empty line
  m_ifs.getline(line, MAX_CHARS); // first receiver-block line; validate
  if (line[0] != '3' || line[1] != '.')
    return -1;
  if (line[2] == 'x')
    return 1; // line 3.x; exit

  if (strncmp("Receiver Type            :", line + 5, 26))
    return -1;
  receiver_type.assign(line + 31, 20);
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Satellite System         :", line + 5, 26))
    return -2;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Serial Number            :", line + 5, 26))
    return -3;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Firmware Version         :", line + 5, 26))
    return -4;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Elevation Cutoff Setting :", line + 5, 26))
    return -5;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Date Installed           :", line + 5, 26)) {
    return -6;
  } else {
    if (str_is_empty(line + 31)) {
      std::cerr << "\n[DEBUG] read_receiver_block: \'Date Installed\' record "
                   "is empty; setting to min date.";
      start = ngpt::datetime<ngpt::milliseconds>::min();
    } else {
      int pos = 31;
      // go to the start of the datetime string
      while (line[pos] && line[pos] == ' ')
        ++pos;
      if (!strncmp(line + pos, "(CCYY-MM-DDThh:mmZ)", 19) ||
          !strncmp(line + pos, "CCYY-MM-DDThh:mmZ", 17)) {
        start = ngpt::datetime<ngpt::milliseconds>::min();
      } else {
        start = strptime_log(line + 31);
      }
    }
  }
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Date Removed             :", line + 5, 26)) {
    return -7;
  } else {
    if (str_is_empty(line + 31)) {
      std::cerr << "\n[DEBUG] read_receiver_block: \'Date Removed\' record "
                   "is empty; setting to max date.";
      stop = ngpt::datetime<ngpt::milliseconds>::max();
    } else {
      int pos = 31;
      // go to the start of the datetime string
      while (line[pos] && line[pos] == ' ')
        ++pos;
      if (!strncmp(line + pos, "(CCYY-MM-DDThh:mmZ)", 19) ||
          !strncmp(line + pos, "CCYY-MM-DDThh:mmZ", 17)) {
        stop = ngpt::datetime<ngpt::milliseconds>::max();
      } else {
        stop = strptime_log(line + 31);
      }
    }
  }
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Temperature Stabiliz.    :", line + 5, 26))
    return -8;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Additional Information   :", line + 5, 26))
    return -9;
  // read additional information lines
  bool new_line = true;
  auto pos = m_ifs.tellg();
  std::size_t nl = 0;
  while (new_line && m_ifs.getline(line, MAX_CHARS)) {
    if (!str_is_empty(line)) {
      new_line = true;
      pos = m_ifs.tellg();
    } else {
      new_line = false;
    }
    if (++nl >= MAX_LINES) {
      return -20;
    }
  }
  m_ifs.seekg(pos);

  return 0;
}

/// Resolve an igs station log's file antenna block (aka the
/// '4.   GNSS Antenna Information' block). This function expects that the
/// instance's stream is ready to read a line of type:
/// '4.1  Antenna Type             : (A20, from rcvr_ant.tab; see
/// instructions)' If the first line read is something different, an error is
/// signaled. It then will read all recorded information and return the
/// important fields. The block should be formatted exactly as described at
/// the IGS specifications (see
/// ftp://igs.org/pub/station/general/sitelog_instr.txt). The function (if
/// successeful), will leave the stream after having read all info of the
/// block, including any 'Additional Information'. That means that the next
/// line to be read, should be an empty line.
///
/// @param[out] antenna_type  The antenna type specified in this block (line
///                           'Antenna Type')
/// @param[out] start         Start time of this block's validity (line
///                           'Date Installed'). If this is set to the string
///                           'CCYY-MM-DDThh:mmZ' inside the block, then
///                           start is set to the min
///                           datetime<ngpt::milliseconds>. The same goes for an
///                           empty string.
/// @param[out] stop          Stop time of this block's validity (line
///                           'Date Removed). If this is set to the string
///                           'CCYY-MM-DDThh:mmZ' inside the block, then
///                           stop is set to the max
///                           datetime<ngpt::milliseconds>. The same goes for an
///                           empty string.
/// @return  The function returns an integer, signaling the following:
///           - int < 0 error
///           - int = 0 all ok
///           - int > 0 block 4.x (skipped; no more blocks to read)
int ngpt::igs_log::read_antenna_block(
    std::string &antenna_type, ngpt::datetime<ngpt::milliseconds> &start,
    ngpt::datetime<ngpt::milliseconds> &stop) {
  constexpr std::size_t MAX_CHARS{256}, MAX_LINES{1000};
  char line[MAX_CHARS];

  m_ifs.getline(line, MAX_CHARS); // empty line
  m_ifs.getline(line, MAX_CHARS); // first antenna-block line; validate
  if (line[0] != '4' || line[1] != '.')
    return -1;
  if (line[2] == 'x')
    return 1; // line 4.x; exit

  if (strncmp("Antenna Type             :", line + 5, 26))
    return -1;
  antenna_type.assign(line + 31, 20);
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Serial Number            :", line + 5, 26))
    return -2;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Antenna Reference Point  :", line + 5, 26))
    return -3;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Marker->ARP Up Ecc. (m)  :", line + 5, 26))
    return -4;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Marker->ARP North Ecc(m) :", line + 5, 26))
    return -5;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Marker->ARP East Ecc(m)  :", line + 5, 26))
    return -6;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Alignment from True N    :", line + 5, 26))
    return -7;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Antenna Radome Type      :", line + 5, 26))
    return -8;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Radome Serial Number     :", line + 5, 26))
    return -9;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Antenna Cable Type       :", line + 5, 26))
    return -10;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Antenna Cable Length     :", line + 5, 26))
    return -11;
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Date Installed           :", line + 5, 26)) {
    return -12;
  } else {
    if (str_is_empty(line + 31)) {
      std::cerr << "\n[DEBUG] read_antenna_block: \'Date Installed\' record "
                   "is empty; setting to min date.";
      start = ngpt::datetime<ngpt::milliseconds>::min();
    } else {
      int pos = 31;
      // go to the start of the datetime string
      while (line[pos] && line[pos] == ' ')
        ++pos;
      if (!strncmp(line + pos, "(CCYY-MM-DDThh:mmZ)", 19) ||
          !strncmp(line + pos, "CCYY-MM-DDThh:mmZ", 17)) {
        start = ngpt::datetime<ngpt::milliseconds>::min();
      } else {
        start = strptime_log(line + 31);
      }
    }
  }
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Date Removed             :", line + 5, 26)) {
    return -13;
  } else {
    if (str_is_empty(line + 31)) {
      std::cerr << "\n[DEBUG] read_antenna_block: \'Date Removed\' record is "
                   "empty; setting to max date.";
      stop = ngpt::datetime<ngpt::milliseconds>::max();
    } else {
      int pos = 31;
      // go to the start of the datetime string
      while (line[pos] && line[pos] == ' ')
        ++pos;
      if (!strncmp(line + pos, "(CCYY-MM-DDThh:mmZ)", 19) ||
          !strncmp(line + pos, "CCYY-MM-DDThh:mmZ", 17)) {
        stop = ngpt::datetime<ngpt::milliseconds>::max();
      } else {
        stop = strptime_log(line + 31);
      }
    }
  }
  if (!m_ifs.getline(line, MAX_CHARS) ||
      strncmp("Additional Information   :", line + 5, 26))
    return -14;
  // read additional information lines
  bool new_line = true;
  auto pos = m_ifs.tellg();
  std::size_t nl = 0;
  while (new_line && m_ifs.getline(line, MAX_CHARS)) {
    if (!str_is_empty(line)) {
      new_line = true;
      pos = m_ifs.tellg();
    } else {
      new_line = false;
    }
    if (++nl >= MAX_LINES) {
      return -20;
    }
  }
  m_ifs.seekg(pos);

  return 0;
}

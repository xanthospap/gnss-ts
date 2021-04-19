#ifndef __NGPT_NOA_EARTHQUAKE_CATALOGUE2_HPP__
#define __NGPT_NOA_EARTHQUAKE_CATALOGUE2_HPP__

///
/// @brief This file defines classes and functions to treat earthquake events
///

#include "earthquake.hpp"
#include "ggdatetime/datetime_read.hpp"
#include "ggdatetime/datetime_write.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include "ggdatetime/dtfund.hpp"
#include "ggeodesy/ellipsoid.hpp"
#include "ggeodesy/geodesy.hpp"
#include "ggeodesy/units.hpp"
#include "ggeodesy/vincenty.hpp"
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace ngpt {

namespace earthquake_catalogue_detail {
/// Max line length (in chars) for an earthquake catalogue file.
constexpr std::size_t MAX_CHARS_IN_LINE = 256;
} // namespace earthquake_catalogue_detail

/// @class earthquake_catalogue
///
/// @brief A wrapper class to hold an earthquake catalogue as distributed by noa
///
/// The use of this class is limited to reading catalogue files; so what it
/// is destined for, is making it easy to read in an earthquake catalogue and
/// producing a list of earthquakes.
///
class earthquake_catalogue {
public:
  /// @brief Constructor using a filename.
  ///
  /// The constructor will try to:
  /// - open the file (using the filename provided) and
  /// - read the file's header and set the stream position (m_end_of_header)
  /// Hence, at exit, the file will be in a read-ready state.
  explicit earthquake_catalogue(const std::string &filename);

  /// No copy constructor.
  earthquake_catalogue(const earthquake_catalogue &) = delete;

  /// No assignment operator.
  earthquake_catalogue &operator=(const earthquake_catalogue &) = delete;

  /// Move constructor.
  earthquake_catalogue(earthquake_catalogue &&ec)
      : m_filename{std::move(ec.m_filename)}, m_ifs{std::move(ec.m_ifs)},
        m_end_of_header{std::move(ec.m_end_of_header)} {}

  /// Move assignment operator.
  earthquake_catalogue &operator=(earthquake_catalogue &&ec) {
    if (this != &ec) {
      m_filename = std::move(ec.m_filename);
      m_ifs = std::move(ec.m_ifs);
      m_end_of_header = std::move(ec.m_end_of_header);
    }
    return *this;
  }

  /// Destructor. If file is open, close it.
  ~earthquake_catalogue() noexcept {
    if (m_ifs.is_open())
      m_ifs.close();
  }

  /// @brief Rewing the stream to the m_end_of_header position.
  ///
  /// Go to the top of the file (just after the header lines); that is, next
  /// line to be read should be the first earthquake in the catalogue.
  void rewind() noexcept { m_ifs.seekg(m_end_of_header, std::ios::beg); }

  /// @brief Read and return the next earthquake
  ///
  /// The function will (try to) read the next earthquake record off from the
  /// catalogue file. If successeful, the eq parameter will hold the (new)
  /// earthquake read, and the function will return true.
  /// If the function fails, then false is returned and eq paramter is not
  /// changed. For the function to fail, it means that EOF is encounterd.
  /// If however the function fails for any other reason (e.g. bad stream,
  /// bad format, etc), then an exception is thrown.
  ///
  /// @param[out] eq  If the function ended successefuly, then eq holds the
  ///                 earthquake read from the input line.
  /// @return         0 if the (next) line was read successefuly and the
  ///                 earthquake was resolved. +1 if EOF was encountered.
  ///                 In case something went wrong -1 is returned.
  int read_next_earthquake(ngpt::earthquake &eq) noexcept;

  /// Navigate to a position in the catalogue file, such that the next
  /// earthquake to be read in will be at an epoch >= to the given epoch
  ///
  /// @param[in] eph  Navigate to the position in file such that the next eq
  ///                 to be read in, occurs at an epoch >= eph
  /// @return The file will set the ifstream buffer to the correct position
  ///         (if found) and return false. +1 will be returned if we reach EOF
  ///         before finding a date >= eph. In case of error, -1 is returned.
  int goto_epoch(ngpt::datetime<ngpt::milliseconds> eph) noexcept;

  void clear_stream_state() { m_ifs.clear(); }

private:
  /// The name of the catalogue
  std::string m_filename;
  /// The input file stream
  std::ifstream m_ifs;
  /// The position within the stream, to start reading the first earthquake
  std::ifstream::pos_type m_end_of_header;

  /// @brief Read header records.
  ///
  /// Header records are the two first lines, which should follow the format:
  /// '      DATE         TIME     LAT.   LONG.  DEPTH    MAGNITUDE'
  /// '                   (GMT)    (N)    (E)    (km)       (Local)'
  /// No line should contain more than
  /// earthquake_catalogue_detail::MAX_CHARS_IN_LINE characters. The function
  /// will also set the m_end_of_header stream position to the point where the
  /// earthquake records start.
  ///
  /// @return True in case the header was read correctly; false otherwise.
  bool read_header() noexcept;

}; // earthquake_catalogue

} // namespace ngpt

#endif

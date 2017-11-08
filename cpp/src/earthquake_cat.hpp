#ifndef __NGPT_EARTHQUAKE_CATALOGUE_HPP__
#define __NGPT_EARTHQUAKE_CATALOGUE_HPP__

///
/// @brief This file defines classes and functions to treat earthquake events
///

// standard headers
#include <fstream>
#include <cstring>
#include <stdexcept>

// datetime headers
#include "ggdatetime/dtcalendar.hpp"
#include "ggdatetime/datetime_read.hpp"

// ggeodesy headers
#include "ggeodesy/ellipsoid.hpp"
#include "ggeodesy/geodesy.hpp"
#include "ggeodesy/vincenty.hpp"

namespace ngpt
{

namespace earthquake_catalogue_detail
{
    /// Max line length (in chars) for an earthquake catalogue file.
    constexpr std::size_t MAX_CHARS_IN_LINE = 256;
}

/// @brief A simple class to hold an earthquake event.
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) will have a time-stamp of type ngpt::datetime<T>.
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    struct earthquake
{
    /// Defult constructor
    earthquake() noexcept
    : m_epoch{},
      m_lon{0},
      m_lat{0},
      m_depth{0},
      m_magnitude{0}
    {}

    /// Full-fledged constructor
    /// @param[in] d    The datetime of the earthquake (when it happened)
    /// @param[in] lat  The latitude of the source (radians)
    /// @param[in] lon  The longtitude of the source (radians)
    /// @param[in] dpth The depth (meters)
    /// @param[in] mag  The magnitude
    explicit
    earthquake(const ngpt::datetime<T>& d, double lat, double lon,
        double dpth, double mag)
    noexcept
    : m_epoch{d},
      m_lon{lon},
      m_lat{lat},
      m_depth{dpth},
      m_magnitude{mag}
    {}

    /// @brief The earthquake's epoch (datetime)
    /// @return A reference to the instance's epoch
    ngpt::datetime<T>& epoch() noexcept { return m_epoch; }

    /// @brief The earthquake's epoch (datetime)
    /// @return The instance's epoch (const version)
    ngpt::datetime<T> epoch() const noexcept { return m_epoch; }

    /// @brief The earthquake's epicenter longtitude.
    /// @return A reference to the instance's longtitude (radians).
    double& lon() noexcept { return m_lon; }
    
    /// @brief The earthquake's epicenter longtitude.
    /// @return The instance's longtitude (radians).
    double lon() const noexcept { return m_lon; }
    
    /// @brief The earthquake's epicenter latitude.
    /// @return A reference to the instance's latitude (radians).
    double& lat() noexcept { return m_lat; }
    
    /// @brief The earthquake's epicenter latitude..
    /// @return The instance's latitude (radians).
    double lat() const noexcept { return m_lat; }
    
    /// @brief The earthquake's magnitude.
    /// @return A reference to the instance's magnitude.
    double& magnitude() noexcept { return m_magnitude; }

    /// @brief The earthquake's magnitude.
    /// @return The instance's magnitude.
    double magnitude() const noexcept { return m_magnitude; }

    /// @brief The earthquake's depth.
    /// @return A reference to the instance's depth (meters).
    double& depth() noexcept { return m_depth; }
    
    /// @brief The earthquake's depth.
    /// @return The instance's depth (meters).
    double depth() const noexcept { return m_depth; }

    /// @brief Return the distance of a point on the ellipsoid from the epcenter.
    ///
    /// The distance is computed along the geodesic that connects the two points
    /// on the ellipsoid E. The line is from the epicenter to the given point.
    /// To compute the distance, the function uses the (inverse) Vincenty
    /// algorithm; hence, the forward and backward azimouths are also computed
    /// and returned (as parameters frw_az and bkw_az).
    ///
    /// @tparam     E       The reference ellipsoid (default is
    ///                     ngpt::ellipsoid::wgs84)
    /// @param[in]  lat     The latitide of the point P (radians)
    /// @param[in]  lon     The longtitude of the point P (radians)
    /// @param[out] frw_az  Forward azimouth (i.e. epicenter to P)
    /// @param[out] bkw_az  Backward azimouth (i.e. P to epicenter)
    /// @return             The distance from the epicenter to point P on the
    ///                     ellipsoid.
    template<ellipsoid E = ellipsoid::wgs84>
        double
        epicenter_distance(double lat, double lon, double& frw_az, double& bkw_az)
        const
    {
        return inverse_vincenty<E>(m_lat, m_lon, lat, lon, frw_az, bkw_az, 1e-12);
    }

private:
    ngpt::datetime<T> m_epoch; ///< The datetime it happened
    double m_lon;              ///< Epicenter longtitude (radians)
    double m_lat;              ///< Epicenter latitude (radians)
    double m_depth;            ///< Depth (meters)
    double m_magnitude;        ///< The magnitude in (??)

}; // end class earthquake

/// A wrapper class to hold an earthquake catalogue as distributed by noa
///
/// The use of this class is limited to reading catalogue files; so what it
/// is destined for, is making it easy to read in an earthquake catalogue and
/// producing a list of earthquakes.
///
/// @tparam T The time precision; this can be any class (of ggdatetime), for
///           which is_of_sec_type is true. This means that T could be e.g.
///           ngpt::seconds, ngpt::milliseconds, etc. The individual epochs
///           (time points) will have a time-stamp of type ngpt::datetime<T>.
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class earthquake_catalogue
{
public:

    /// @brief Constructor using a filename.
    ///
    /// The constructor will try to:
    /// - open the file (using the filename provided) and
    /// - read the file's header and set the stream position (m_end_of_header)
    /// Hence, at exit, the file will be in a read-ready state.
    explicit earthquake_catalogue(const std::string& filename)
    : m_filename{filename},
      m_ifs{filename.c_str(), std::ifstream::in},
      m_end_of_header{0}
    {
        if ( !m_ifs.is_open() ) {
            throw std::invalid_argument
            ("earthquake_catalogue: Could not open earthquake catalogue file: \""+filename+"\"");
        }
        if ( !read_header() ) {
            m_ifs.close();
            throw std::invalid_argument
            ("earthquake_catalogue: Could not read header in earthquake catalogue file: \""+filename+"\"");
        }

    }
    
    /// No copy constructor.
    earthquake_catalogue(const earthquake_catalogue&) = delete;

    /// No assignment operator.
    earthquake_catalogue& operator=(const earthquake_catalogue&) = delete;

    /// Move constructor.
    earthquake_catalogue(earthquake_catalogue&& ec)
    : m_filename{std::move(ec.m_filename)},
      m_ifs{std::move(ec.m_ifs)},
      m_end_of_header{std::move(ec.m_end_of_header)}
    {}

    /// Move assignment operator.
    earthquake_catalogue& operator=(earthquake_catalogue&& ec)
    {
        if (*this != ec) {
            m_filename = std::move(ec.m_filename);
            m_ifs = std::move(ec.m_ifs);
            m_end_of_header = std::move(ec.m_end_of_header);
        }
        return *this;
    }

    /// Destructor. If file is open, close it.
    ~earthquake_catalogue() noexcept
    { 
        if (m_ifs.is_open()) { m_ifs.close(); }
    }

    /// @brief Rewing the stream to the m_end_of_header position.
    ///
    /// Go to the top of the file (just after the header lines); that is, next
    /// line to be read should be the first earthquake in the catalogue.
    void rewind() noexcept
    {
        m_ifs.seekg(m_end_of_header, std::ios::beg);
        return;
    }

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
    /// @return         True if the (next) line was read successefuly and the
    ///                 earthquake was resolved. False if EOF was encountered.
    ///                 In case something went wrong, an exception is thrown.
    bool read_next_earthquake(earthquake<T>& eq)
    {
        static char line[earthquake_catalogue_detail::MAX_CHARS_IN_LINE];

        if ( !m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE) )
        {
            if ( !m_ifs.eof() ) {
                throw std::invalid_argument(
                "read_next_earthquake: Error reading earthquake catalogue file: \""+m_filename+"\".");
            }
            return false;
        }
        
        // constexpr double deg2rad = ngpt::DPI / 180.0;
        static float info[4];
        char *start(line), *end;
        datetime<T> eph = ngpt::strptime_yod_hms<T>(line, &start);

        for (int i = 0; i < 4; ++i) {
            info[i] = std::strtod(start, &end);
            if (errno == ERANGE || start == end) {
                errno = 0;
                throw std::invalid_argument
                    ("read_next_earthquake: Invalid line: \""+std::string(line)+"\" (argument #"+
                    std::to_string(i+1)+") in catalogue file.");
            }
            start = end;
        }

        earthquake<T> eqt {eph, deg2rad(info[0]), deg2rad(info[1]), info[2]/1000.0, info[3]};
        eq = std::move(eqt);
        
        return true;
    }
    
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
    /// No line should contain more than earthquake_catalogue_detail::MAX_CHARS_IN_LINE
    /// characters. The function will also set the m_end_of_header stream
    /// position to the point where the earthquake records start.
    ///
    /// @return True in case the header was read correctly; false otherwise.
    bool read_header() noexcept
    {
        char line[earthquake_catalogue_detail::MAX_CHARS_IN_LINE];
        const char* line1 = "      DATE         TIME     LAT.   LONG.  DEPTH    MAGNITUDE";
        const char* line2 = "                   (GMT)    (N)    (E)    (km)       (Local)";

        if ( !m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE)
            || std::strncmp(line, line1, std::strlen(line1)) ) 
        {
            return false;
        }
        if ( !m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE)
            || std::strncmp(line, line2, std::strlen(line2)) )
        {
            return false;
        }
        m_end_of_header = m_ifs.tellg();
        return true;
    }

}; // end class earthquake_catalogue

} // end namespace ngpt

#endif

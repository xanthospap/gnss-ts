#ifndef __EARTHQUAKE_CATALOGUE__
#define __EARTHQUAKE_CATALOGUE__

#include <fstream>
#include <cstring>
#include <stdexcept>
#include "dtcalendar.hpp"
#include "datetime_read.hpp"
#include "geoconst.hpp"   // for DPI

namespace ngpt
{

/// Max line length for catalogue files
namespace earthquake_catalogue_detail
{
    constexpr std::size_t MAX_CHARS_IN_LINE = 256;
}

/// A simple class to hold an earthquake info.
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    struct earthquake
{
    /// Empty!
    earthquake() noexcept
    : m_date{}, m_longtitude{0}, m_latitude{0}, m_depth{0}, m_magnitude{0}
    {}

    /// Full-fledged constructor
    explicit earthquake(const ngpt::datetime<T>& d, double lon, double lat,
        double dpth, double mag) noexcept
    : epoch{d}, longtitude{lon}, latitude{lat}, depth{dpth},
      magnitude{mag}
    {}

    /// the date it appeared
    ngpt::datetime<T> epoch;
    
    /// The longtitude, latitude and depth (in radians, radians, meters)
    double longtitude, latitude, depth;
    
    /// The magnitude in (??)
    double magnitude;

}; // end class earthquake

/// A wrapper class to hold an earthquake catalogue as distributed by noa
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    class earthquake_catalogue
{
public:

    /// Constructor using a filename.
    earthquake_catalogue(const std::string& filename)
    : m_filename{filename},
      m_ifs{filename.c_str(), std::ifstream::in},
      m_end_of_header{0}
    {
        if ( !m_ifs.is_open() ) {
            throw std::invalid_argument
            ("Could not open earthquake catalogue file: \""+filename+"\"");
        }
        if ( !read_header() ) {
            m_ifs.close();
            throw std::invalid_argument
            ("Could not read header in earthquake catalogue file: \""+filename+"\"");
        }

    }
    
    /// No copy constructor
    earthquake_catalogue(const earthquake_catalogue&) = delete;

    /// No assignment operator
    earthquake_catalogue& operator=(const earthquake_catalogue&) = delete;

    /// Move constructor
    earthquake_catalogue(earthquake_catalogue&& ec)
    : m_filename{std::move(ec.m_filename)},
      m_ifs{std::move(ec.m_ifs)},
      m_end_of_header{std::move(ec.m_end_of_header)}
    {}

    /// Move assignment operator
    earthquake_catalogue& operator=(earthquake_catalogue&& ec)
    {
        if (*this != ec) {
            m_filename = std::move(ec.m_filename);
            m_ifs = std::move(ec.m_ifs);
            m_end_of_header = std::move(ec.m_end_of_header);
        }
        return *this;
    }

    /// Destructor
    ~earthquake_catalogue() noexcept
    { 
        if (m_ifs.is_open()) { m_ifs.close(); }
    }

    /// Go to the top of the file (after the header lines); that is, next line
    /// to be read should be the first earthquake in the catalogue.
    void rewind() noexcept
    {
        m_ifs.seekg(m_end_of_header, std::ios::beg);
        return;
    }

    /// Read and return the next earthquake
    bool read_next_earthquake(earthquake<T>& eq)
    {
        static char line[earthquake_catalogue_detail::MAX_CHARS_IN_LINE];

        if ( !m_ifs.getline(line, earthquake_catalogue_detail::MAX_CHARS_IN_LINE) ) 
            return false;
        
        constexpr double deg2rad = ngpt::DPI / 180.0;
        static float info[4];
        char *start(line), *end;
        datetime<T> eph =  ngpt::strptime_yod_hms<T>(line, &start);
        // std::cout<<"\tAfter reading date, start is at: "<< start - line << " dld: " << start << "\n";

        for (int i = 0; i < 4; ++i) {
            info[i] = std::strtod(start, &end);
            if (errno == ERANGE || start == end) {
                errno = 0;
                throw std::invalid_argument
                    ("Invalid line: \""+std::string(line)+"\" (argument #"+
                    std::to_string(i+1)+") in catalogue file.");
            }
            start = end;
        }
        earthquake<T> eqt {eph, info[0]*deg2rad, info[1]*deg2rad, info[2]/1000.0, info[3]};
        eq = eqt;
        return true;
    }
    
private:
    /// the name of the catalogue
    std::string m_filename;
    /// the input file stream
    std::ifstream m_ifs;
    /// the position within the stream, to start reading the first earthquake
    std::ifstream::pos_type m_end_of_header;
    
    /// read header records
    bool read_header()
    {
        char line[256];
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

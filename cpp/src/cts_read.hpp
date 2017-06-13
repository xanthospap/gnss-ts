#ifndef __NGPT_CRD_TSREAD__
#define __NGPT_CRD_TSREAD__

// standard headers
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cctype>

// ggdatetime headers
#include "ggdatetime/datetime_read.hpp"
#include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "crdts.hpp"
// #include "genflags.hpp"
#include "tsflag.hpp"
#include "crdts.hpp"

namespace ngpt
{

/// Conver a c-string to a flag<pt_marker>. Leading whitespace characters are
/// ignored, untill a non-whitespace char is encountered.
flag<pt_marker>
str2flag(const char* str, char** stop)
{
    const char* c = str;
    flag<pt_marker> f;

    while ( *c && (*c == ' ') ) ++c; // go to the first non-whitespace char

    while ( *c && (*c != ' ') ) {
        if ( std::isdigit(*c) || *c == '-' ) return f;
        if (*c == 's') f.set(pt_marker::skip);
        else if (*c == 'o') f.set(pt_marker::outlier);
        else throw 1;// not very informative .... TODO
        ++c;
    }
    *stop = (char*)c; // fuck it, just cast to non-const
    return f;
}

/// Read a time-series file, as fromated by the function:
/// crdts<T>.dump(...)
///
template<typename T,
         typename = std::enable_if_t<T::is_of_sec_type>
         >
    crdts<T> readin(const std::string& infile, const std::string& ts_name)
{
    constexpr std::size_t MAX_CHARS {256};

    std::ifstream ifs (infile.c_str(), std::ifstream::in);
    if ( !ifs.is_open() ) {
        throw std::invalid_argument("Cannot find file \""+infile+"\"");
    }

    char line[MAX_CHARS];
    char *cptr(&line[0]), *end;
    double data[7], mjd, fmjd;
    flag<pt_marker> fx, fy, fz;
    crdts<T> ts {ts_name};
#ifdef DEBUG
    std::size_t line_counter = 0;
#endif

    while ( ifs.getline(line, MAX_CHARS) ) {
        if ( *line != '#' ) {
            cptr = line;
            for (int i = 0; i < 3; ++i) {
                data[i] = std::strtod(cptr, &end);
                if (errno == ERANGE) {
                    errno = 0;
                    ifs.close();
                    throw std::invalid_argument
                        ("Invalid record line: \""+std::string(line)+"\" (argument #"+std::to_string(i)+")");
                }
                cptr = end;
            }
            fx = str2flag(cptr, &end);
            cptr = end;

            for (int i = 3; i < 5; ++i) {
                data[i] = std::strtod(cptr, &end);
                if (errno == ERANGE) {
                    errno = 0;
                    ifs.close();
                    throw std::invalid_argument
                        ("Invalid record line: \""+std::string(line)+"\" (argument #"+std::to_string(i)+")");
                }
                cptr = end;
            }
            fy = str2flag(cptr, &end);
            cptr = end;

            for (int i = 5; i < 7; ++i) {
                data[i] = std::strtod(cptr, &end);
                if (errno == ERANGE) {
                    errno = 0;
                    ifs.close();
                    throw std::invalid_argument
                        ("Invalid record line: \""+std::string(line)+"\" (argument #"+std::to_string(i)+")");
                }
                cptr = end;
            }
            fz = str2flag(cptr, &end);
            cptr = end;
            
            fmjd = std::modf(data[0], &mjd);
            fmjd *= static_cast<double>(T::max_in_day);
            datetime<T> epoch {modified_julian_day{static_cast<long>(mjd)},
                T{static_cast<long>(fmjd)}};
            ts.add(epoch, data[1], data[3], data[5], data[2], data[4], data[6], fx, fy, fz);
#ifdef DEBUG
            ++line_counter;
#endif
        }
    }

#ifdef DEBUG
    std::cout<<"\tRead #"<<line_counter<<" lines from cts file.\n";
#endif
    ifs.close();
    ts.crd_type() = coordinate_type::unknown;
    return ts;
}

/// Read a time-series in cartesian coordinates (aka crdts<T>) off from a .cts file.
/// Note that each line (in the file) can have no more than 256 characters.
/// Also note that the sigmas are scaled to 1000 (i.e. they are assumed meters
/// and converted to millimeters).
template<class T,
        typename = std::enable_if_t<T::is_of_sec_type>
        >
    crdts<T> cts_read(const std::string& cts_file, const std::string& ts_name)
{
    constexpr std::size_t MAX_CHARS {256};

    std::ifstream ifs (cts_file.c_str(), std::ifstream::in);
    if ( !ifs.is_open() ) {
        throw std::invalid_argument("Cannot find file \""+cts_file+"\"");
    }

    char line[MAX_CHARS];
    char *cptr(&line[0]), *end;
    ngpt::datetime<T> epoch;
    double data[6];
    crdts<T> ts {ts_name};
#ifdef DEBUG
    std::size_t line_counter = 0;
#endif

    while ( ifs.getline(line, MAX_CHARS) ) {
        if ( *line != '#' ) {
            epoch = ngpt::strptime_ymd_hms<T>(line, &cptr);
            ++cptr;
            for (int i = 0; i < 6; ++i) {
                data[i] = std::strtod(cptr, &end);
                if (errno == ERANGE) {
                    errno = 0;
                    ifs.close();
                    throw std::invalid_argument
                        ("Invalid record line: \""+std::string(line)+"\" (argument #"+std::to_string(i)+")");
                }
                cptr = end;
            }
            ts.add(epoch, data[0], data[2], data[4], data[1]*1000.0,
                data[3]*1000.0, data[5]*1000.0);
#ifdef DEBUG
            ++line_counter;
#endif
        }
    }

#ifdef DEBUG
    std::cout<<"\tRead #"<<line_counter<<" lines from cts file.\n";
#endif
    ifs.close();
    ts.crd_type() = coordinate_type::cartesian;
    return ts;
}

} // end namespace ngpt

#endif

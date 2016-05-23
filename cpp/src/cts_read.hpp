#ifndef __NGPT_CRD_TSREAD__
#define __NGPT_CRD_TSREAD__

#include <iostream>
#include <cstdlib>
#include <fstream>
#include "crdts.hpp"
#include "dtcalendar.hpp"
#include "genflags.hpp"
#include "tsflagenum.hpp"
#include "timeseries.hpp"
#include "crdts.hpp"

namespace ngpt
{

/// Read a time-series (aka crdts<T>) off from a .cts file.
/// Note that each line (in the file) can have no more than 256 characters
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
    char *cptr, *end;
    ngpt::datetime<T> epoch;
    double data[6];

    crdts<T> ts{ts_name};

    while ( ifs.getline(line, MAX_CHARS) ) {
        if ( *line != '#' ) {
            epoch = ngpt::strptime_ymd_hms(line, cptr);
            ++cptr;
            for (int i = 0; i < 6; ++i) {
                data[i] = std::strtod(cptr, &end);
                if (errno == ERANGE) {
                    errno = 0;
                    ifs.close();
                    throw std::invalid_argument
                        ("Invalid record line: \""+std::string(line)+"\" (argument #" << std::to_string(i)<<")");
                }
            }
            ts.add(epoch, data[0], data[2], data[4], data[1], data[3], data[5]);
        }
    }

    ifs.close();
    return ts;
}

} // end namespace ngpt

#endif

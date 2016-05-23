#ifndef __NGPT_DT_READERS__
#define __NGPT_DT_READERS__

#include <cstdlib>
#include "dtfund.hpp"
#include "dtcalendar.hpp"
#ifdef DEBUG
    #include <iostream>
#endif

namespace ngpt {

/// Time format
enum class datetime_format : char {
    ymd,
    ymd_hms,
    ydoy,
    ydoy_hms
};

/// Resolve a date from a c-string.
/*
template<typename T,
        typename = std::enable_if_t<T::is_of_sec_type>,
        datetime_format F
        >
    datetime<T> strptime(const char* str)
{}
*/

/// Read a date from a c-string of type:
/// YYYY-MM-DD HH:MM:SS.SSSS
template<typename T> datetime<T> strptime_ymd_hms(const char* str)
{
    char *end;
    const char* start = str;
    int ints[5];
    double secs;

    for (int i = 0; i < 5; ++i) {
        ints[i] = static_cast<int>( std::abs(std::strtol(start, &end, 10)) );
        // std::cout << "\tArgument " << i << " = " << ints[i] << ", string: [" << std::string(start, end-start) << "]\n";
        if (errno == ERANGE || start == end) {
            // if (errno == ERANGE) {std::cout << "\tERRNO set!\n";}
            // if (start == end) {std::cout<<"\tstart == end\n";}
            errno = 0;
            throw std::invalid_argument
                ("Invalid date format: \""+std::string(str)+"\" (argument " + std::to_string(i+1) + ").");
        }
        start = end+1;
    }
    secs = std::strtod(start, &end);
    if (errno == ERANGE) {
        errno = 0;
        throw std::invalid_argument
            ("Invalid date format: \""+std::string(str)+"\" (argument 6)");
    }
    return datetime<T> {year{ints[0]}, month{ints[1]}, day_of_month{ints[2]},
        hours{ints[3]}, minutes{ints[4]}, secs};
}

} // namespace ngpt

#endif

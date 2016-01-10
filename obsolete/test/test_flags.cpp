#include <iostream>
#include "flags.hpp"

using ts::ts_flag;

int main ()
{
    // initialize a flag; all fields should be '0'
    ts::flag_type flag;
    std::cout << "\nDefault initialized flag: ";
    ts::bs_print( flag );

    // check if it is an outlier
    if ( check_flag(flag, ts::ts_flag::outlier) ) {
        std::cout << "\nThis is an outlier!";
    } else {
        std::cout << "\nNo, not an outlier!";
    }

    // set this flag to be an outlier and print it
    std::cout << "\nSetting the outlier flag:";
    flag |= ts::ts_flag::outlier;
    ts::bs_print( flag );

    // also set velocity change and earthquake
    std::cout << "\nSetting the vel. change and earthquake flags:";
    flag |= (ts::ts_flag::velchg | ts::ts_flag::ethq);
    ts::bs_print( flag );

    // check multiple flags
    std::cout << "\nChecking multiple flags (outlier, skip, vel. change):";
    if (  check_flag(flag, ts::ts_flag::outlier)
       || check_flag(flag, ts::ts_flag::skip)
       || check_flag(flag, ts::ts_flag::velchg)
       ) {
        std::cout << "\nOne (or more) of them flags set!";
    } else {
        std::cout << "\nNone of them flags set!";
    }

    std::cout << "\n";
    return 0;
}

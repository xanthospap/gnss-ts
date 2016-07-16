#include <iostream>

#include "cts_read.hpp"
#include "crdts.hpp"
#include "genflags.hpp"

//
// This is a test program to check the ngpt::crdts<> class.
//

int
main(int argc, char* argv[])
{
    if (argc < 2 || argc > 3) {
        std::cerr<<"Usage: read_cts <cts file> [catalogue]\n";
        return 1;
    }

    std::string cts_file = std::string(argv[1]);
    std::string cts_name = std::string("test");

    // Read in the coordinate time-series from the input file
    std::cout << "Reading time-series file \"" <<cts_file<<"\", as station \""<<cts_name<<"\"\n";
    ngpt::crdts<ngpt::milliseconds> ts
        = ngpt::cts_read<ngpt::milliseconds>(cts_file, cts_name);

    // Transform cartesian to topocentric (including sigmas)
    ts.cartesian2topocentric();

    // If catalogue file provided, read and apply the interesting earthquakes
    if (argc > 2) {
        ngpt::earthquake_catalogue<ngpt::milliseconds> eq_cat
            {std::string(argv[2])};
        ts.apply_earthquake_catalogue(eq_cat);
    }

    return 0;
}

#include <fstream>
#include "tmsr.hpp"

using ts::time_series;

std::size_t
time_series::read_from_cts(const char* filename)
{
    std::ifstream fin { filename };
    if ( !fin.is_open() ) {
        return 0;
    }

    char line[256];


}

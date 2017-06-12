#include <iostream>
#include "timeseries.hpp"

using ngpt::datetime;
using ngpt::milliseconds;

int
main(int argc, char* argv[])
{
    // Let's construct a vector of epochs, to be approximately of size
    // 15(years) * (365 days/year) = 5475, let's round it to 6000
    int NUM_EPOCHS = 6000;
    long MJD0      = 53743;
    std::vector<datetime<milliseconds>> tvec;
    for (int i = 0; i < NUM_EPOCHS; i++) tvec.emplace_back();

    std::cout<<"\n";
    return 0;
}

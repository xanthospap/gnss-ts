#include <iostream>
#include "stalogrd.hpp"
#include "ggdatetime/dtcalendar.hpp"

using ngpt::seconds;
int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::cerr<<"\nUsage: stalogrd <log_file>\n";
        return 1;
    }

    ngpt::igs_log log { std::string(argv[1]) };
    log.receiver_changes<seconds>();

    return 0;
}

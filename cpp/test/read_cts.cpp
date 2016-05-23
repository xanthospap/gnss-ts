#include <iostream>
#include "cts_read.hpp"

int main(int argc, char* argv[])
{
    if (argc != 1) {
        std::cerr<<"Usage: read_cts <cts file>\n";
        return 1;
    }

    std::string cts_file = std::string(argv[1]);
    std::string cts_name = "test";
    ngpt::crdts<ngpt::milliseconds> ts = ngpt::cts_read<ngpt::milliseconds>(cts_file, cts_name);

    return 0;
}

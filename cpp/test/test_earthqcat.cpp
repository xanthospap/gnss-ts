#include "dtfund.hpp"
#include "earthquake_cat.hpp"

using ngpt::seconds;

int main(int argc, char* argv[])
{

    if (argc != 2) {
        std::cerr<<"Usage: read_cat <cat file>\n";
        return 1;
    }

    std::string cat_file = std::string(argv[1]);
    std::cout << "Reading catalogue file \"" <<cat_file<<"\"n";

    ngpt::earthquake<seconds> eq;
    ngpt::earthquake_catalogue<seconds> ecat {cat_file};
    std::size_t eqk_num = 0;

    while ( ecat.read_next_earthquake(eq) ) {
        ++eqk_num;
    }

    std::cout << "Read #"<<eqk_num<<" earthquakes\n";
    return 0;
}

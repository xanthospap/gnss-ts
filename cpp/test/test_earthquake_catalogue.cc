#include "earthquake_catalogue.hpp"
#include <iostream>
#include <string>

using namespace ngpt;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "usage: $>" << argv[0] << " <noa_earthquake_catalogue_file>\n";
    return 1;
  }

  int status;

  // create an earthquake_catalogue instance
  earthquake_catalogue cat{std::string(argv[1])};

  // an earthquake instance
  earthquake eq;

  // get and print the first earthquake in the file and print it
  assert(!(status = cat.read_next_earthquake(eq)));
  std::cout << "First earthquake in file is: " << eq.to_string() << "\n";

  // random date
  datetime<milliseconds> rand_t(year(2007), month(9), day_of_month(20),
                                hours(23), minutes(1), milliseconds(13000));

  // go to random epoch in file, and print the first earthquake after that
  assert(!cat.goto_epoch(rand_t));
  assert(!cat.read_next_earthquake(eq));
  assert(eq.m_epoch >= rand_t);
  std::cout << "First earthquake after 2007/09/20 23:01 is: " << eq.to_string()
            << "\n";

  // random late date ....
  datetime<milliseconds> rand_t2(year(2071), month(9), day_of_month(20),
                                 hours(23), minutes(1), milliseconds(13000));
  // go to random epoch in file, and print the first earthquake after that
  status = cat.goto_epoch(rand_t2);
  // well.... we should have hit EOF
  assert(status > 0);
  // so that the follwing returns an error
  status = cat.read_next_earthquake(eq);
  assert(status > 0);

  std::cout << "\n";
  return 0;
}

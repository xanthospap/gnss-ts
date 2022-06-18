#include "datetime/datetime_write.hpp"
#include "igs_sta_log.hpp"
#include <iostream>

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "usage: $>" << argv[0] << " <igs_log_file>\n";
    return 1;
  }

  // int status;

  // create an igs_log instance
  igs_log log{std::string(argv[1])};

  // get receiver changes
  auto recvec = log.receiver_changes();
  std::cout << "Receiver Changes: \n";
  for (const auto &d : recvec)
    std::cout << "\t->" << strftime_ymd_hmfs(d) << "\n";

  // get antenna changes
  auto antvec = log.antenna_changes();
  std::cout << "Antenna Changes: \n";
  for (const auto &d : antvec)
    std::cout << "\t->" << strftime_ymd_hmfs(d) << "\n";

  std::cout << "\n";
  return 0;
}

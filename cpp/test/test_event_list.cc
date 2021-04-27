#include "event_list.hpp"
#include <cassert>
#include <iostream>

bool validate_event_list(const ngpt::event_list &e) {
  auto it_prev = e.it_cbegin(), it_next = e.it_cbegin() + 1;
  while (it_next != e.it_cend()) {
    assert(it_prev->epoch() <= it_next->epoch());
    if (it_prev->epoch() == it_next->epoch()) {
      assert(it_prev->event_type() != it_next->event_type());
    }
    ++it_prev;
    ++it_next;
  }
  return true;
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <event list file> <igs log file>\n";
    return 1;
  }

  // create an event_list instance
  ngpt::event_list list;

  // parse and apply events recorded in the event list file
  list.apply_event_list_file(argv[1]);

  // parse and apply events recorded in the station log file
  list.apply_stalog_file(argv[2]);

  // write events to STDOUT
  list.dump_event_list(std::cout);

  // validate instance
  assert(validate_event_list(list));

  std::cout << "\n";
  return 0;
}

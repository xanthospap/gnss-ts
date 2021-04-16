#include "ts_flag.hpp"
#include <cassert>

using ngpt::data_flag;

int main() {

  // construct an instance; should be clean
  data_flag f;
  assert(f.is_clean());

  // set it to an outlier
  f.set(ngpt::pt_marker::outlier);
  assert(f.is_set(ngpt::pt_marker::outlier) && !f.is_set(ngpt::pt_marker::skip) && !f.is_clean());

  // toggle skip flag on
  f.set(ngpt::pt_marker::skip);
  assert(f.is_set(ngpt::pt_marker::outlier) && f.is_set(ngpt::pt_marker::skip) && !f.is_clean());

  // un-set outlier
  f.clear(ngpt::pt_marker::outlier);
  assert(!f.is_set(ngpt::pt_marker::outlier) && f.is_set(ngpt::pt_marker::skip) && !f.is_clean());

  // un-set skip
  f.clear(ngpt::pt_marker::skip);
  assert(!f.is_set(ngpt::pt_marker::outlier) && !f.is_set(ngpt::pt_marker::skip) && f.is_clean());

  // by now we should have un-set everything
  assert(f == data_flag{});
  
  // set all, and clear
  f.clear(ngpt::pt_marker::outlier);
  f.set(ngpt::pt_marker::skip);
  f.clear();
  assert(f == data_flag{});

  return 0;
}

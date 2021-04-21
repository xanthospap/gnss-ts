#include "genflags.hpp"
#include "ts_flag.hpp"
#include <cassert>

int main() {

  using pt_marker_flag = ngpt::flag<ngpt::pt_marker>;

  // construct an instance; should be clean
  // ngpt::flag<ngpt::pt_marker> f;
  pt_marker_flag f;
  assert(f.is_clean());

  // set it to an outlier
  f.set(ngpt::pt_marker::outlier);
  assert(f.is_set(ngpt::pt_marker::outlier) &&
         !f.is_set(ngpt::pt_marker::skip) && !f.is_clean());

  // toggle skip flag on
  f.set(ngpt::pt_marker::skip);
  assert(f.is_set(ngpt::pt_marker::outlier) &&
         f.is_set(ngpt::pt_marker::skip) && !f.is_clean());

  // un-set outlier
  f.clear(ngpt::pt_marker::outlier);
  assert(!f.is_set(ngpt::pt_marker::outlier) &&
         f.is_set(ngpt::pt_marker::skip) && !f.is_clean());

  // un-set skip
  f.clear(ngpt::pt_marker::skip);
  assert(!f.is_set(ngpt::pt_marker::outlier) &&
         !f.is_set(ngpt::pt_marker::skip) && f.is_clean());

  // by now we should have un-set everything
  assert(f == ngpt::flag<ngpt::pt_marker>{});

  // set all, and clear
  f.clear(ngpt::pt_marker::outlier);
  f.set(ngpt::pt_marker::skip);
  f.clear();
  assert(f == ngpt::flag<ngpt::pt_marker>{});

  return 0;
}

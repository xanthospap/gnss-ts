#include "genflags.hpp"
#include "ts_flag.hpp"
#include <cassert>

int main() {

  using pt_marker_flag = dso::flag<dso::pt_marker>;

  // construct an instance; should be clean
  // dso::flag<dso::pt_marker> f;
  pt_marker_flag f;
  assert(f.is_clean());

  // set it to an outlier
  f.set(dso::pt_marker::outlier);
  assert(f.is_set(dso::pt_marker::outlier) &&
         !f.is_set(dso::pt_marker::skip) && !f.is_clean());

  // toggle skip flag on
  f.set(dso::pt_marker::skip);
  assert(f.is_set(dso::pt_marker::outlier) &&
         f.is_set(dso::pt_marker::skip) && !f.is_clean());

  // un-set outlier
  f.clear(dso::pt_marker::outlier);
  assert(!f.is_set(dso::pt_marker::outlier) &&
         f.is_set(dso::pt_marker::skip) && !f.is_clean());

  // un-set skip
  f.clear(dso::pt_marker::skip);
  assert(!f.is_set(dso::pt_marker::outlier) &&
         !f.is_set(dso::pt_marker::skip) && f.is_clean());

  // by now we should have un-set everything
  assert(f == dso::flag<dso::pt_marker>{});

  // set all, and clear
  f.clear(dso::pt_marker::outlier);
  f.set(dso::pt_marker::skip);
  f.clear();
  assert(f == dso::flag<dso::pt_marker>{});

  return 0;
}

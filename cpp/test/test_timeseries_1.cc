#include "timeseries.hpp"
#include <datetime/datetime_write.hpp>
#include <datetime/dtcalendar.hpp>
#include <random>

std::default_random_engine e;
std::uniform_real_distribution<> dis(0, 1);

void dummy_print(const dso::timeseries &ts) {
  std::size_t idx = 0;
  while (idx < ts.size()) {
    std::cout << dso::strftime_ymd_hmfs(ts.epochs_vec()[idx]) << " "
              << ts.data_points_vec()[idx].value() << "\n";
    ++idx;
  }
  return;
}

int main() {

  // vector of epochs
  std::vector<dso::datetime<dso::milliseconds>> epochs;
  dso::datetime<dso::milliseconds> t{dso::year(2019), dso::day_of_year(1)};
  dso::datetime_interval<dso::milliseconds> dt{dso::modified_julian_day(1),
                                                 dso::milliseconds(1000)};
  while (t < dso::datetime<dso::milliseconds>{dso::year(2020),
                                                dso::day_of_year(1)}) {
    epochs.emplace_back(t);
    t += dt;
  }

  // create a time-series instance
  dso::timeseries ts(&epochs);
  std::size_t idx = 0;
  while (idx++ < epochs.size())
    ts.add_point(dis(e));

  // copy the timeseries
  {
    dso::timeseries ts2(ts);
    assert(ts2 == ts);
    ts2.data_points_vec()[1] = 999.99e0;
    assert(ts2 != ts);
  }

  // let's mark as outliers everything in .... February 2019
  auto t1 = dso::datetime<dso::milliseconds>{dso::year(2019), dso::month(2),
                                               dso::day_of_month(1)};
  auto t2 = dso::datetime<dso::milliseconds>{dso::year(2019), dso::month(3),
                                               dso::day_of_month(1)};
  auto marked_pts = ts.mark(dso::pt_marker::outlier, t1, t2);
  assert(marked_pts == 28);
  idx = 0;
  while (idx < epochs.size()) {
    t = ts.epochs_vec()[idx];
    const auto p = ts.data_points_vec()[idx];
    if (t < t1 || t >= t2) {
      assert(p.flag().is_clean());
    } else {
      assert(p.flag().is_set(dso::pt_marker::outlier));
    }
    ++idx;
  }

  // let's cut the time series ....
  auto new_epochs = epochs;
  dso::timeseries ts2(ts);
  auto new_size = ts2.cut(t1, t2, &new_epochs);
  ts2.set_epoch_vec(&new_epochs);
  assert((std::size_t)new_size == new_epochs.size() &&
         new_epochs.size() == ts2.data_points_vec().size());
  assert(new_epochs[0] >= t1 && new_epochs[new_size - 1] < t2);
  for (int i = 1; i < 12; i++) {
    auto dt1 = dso::datetime<dso::milliseconds>{
        dso::year(2019), dso::month(i), dso::day_of_month(1)};
    auto dt2 = dso::datetime<dso::milliseconds>{
        dso::year(2019), dso::month(i + 1), dso::day_of_month(1)};
    new_epochs = epochs;
    ts2 = ts;
    new_size = ts2.cut(dt1, dt2, &new_epochs);
    ts2.set_epoch_vec(&new_epochs);
    assert((std::size_t)new_size == new_epochs.size() &&
           new_epochs.size() == ts2.data_points_vec().size());
    assert(new_epochs[0] >= dt1 && new_epochs[new_size - 1] < dt2);
    assert(ts2.data_points_vec().size() == ts2.epochs_vec().size());
  }
  // dummy_print(ts2);

  // print it
  // dummy_print(ts);

  std::cout << "\n";
  return 0;
}

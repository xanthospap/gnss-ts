#include "timeseries.hpp"
#include <ggdatetime/datetime_write.hpp>
#include <ggdatetime/dtcalendar.hpp>
#include <random>

std::default_random_engine e;
std::uniform_real_distribution<> dis(0, 1);

void dummy_print(const ngpt::timeseries &ts) {
  std::size_t idx = 0;
  while (idx < ts.size()) {
    std::cout << ngpt::strftime_ymd_hmfs(ts.epochs_vec()[idx]) << " "
              << ts.data_points_vec()[idx].value() << "\n";
    ++idx;
  }
  return;
}

int main() {

  // vector of epochs
  std::vector<ngpt::datetime<ngpt::milliseconds>> epochs;
  ngpt::datetime<ngpt::milliseconds> t{ngpt::year(2019), ngpt::day_of_year(1)};
  ngpt::datetime_interval<ngpt::milliseconds> dt{ngpt::modified_julian_day(1),
                                                 ngpt::milliseconds(1000)};
  while (t < ngpt::datetime<ngpt::milliseconds>{ngpt::year(2020),
                                                ngpt::day_of_year(1)}) {
    epochs.emplace_back(t);
    t += dt;
  }

  // create a time-series instance
  ngpt::timeseries ts(&epochs);
  std::size_t idx = 0;
  while (idx++ < epochs.size())
    ts.add_point(dis(e));

  // copy the timeseries
  {
    ngpt::timeseries ts2(ts);
    assert(ts2 == ts);
    ts2.data_points_vec()[1] = 999.99e0;
    assert(ts2 != ts);
  }

  // let's mark as outliers everything in .... February 2019
  auto t1 = ngpt::datetime<ngpt::milliseconds>{ngpt::year(2019), ngpt::month(2),
                                               ngpt::day_of_month(1)};
  auto t2 = ngpt::datetime<ngpt::milliseconds>{ngpt::year(2019), ngpt::month(3),
                                               ngpt::day_of_month(1)};
  auto marked_pts = ts.mark(ngpt::pt_marker::outlier, t1, t2);
  assert(marked_pts == 28);
  idx = 0;
  while (idx < epochs.size()) {
    t = ts.epochs_vec()[idx];
    const auto p = ts.data_points_vec()[idx];
    if (t < t1 || t >= t2) {
      assert(p.flag().is_clean());
    } else {
      assert(p.flag().is_set(ngpt::pt_marker::outlier));
    }
    ++idx;
  }

  // let's cut the time series ....
  auto new_epochs = epochs;
  ngpt::timeseries ts2(ts);
  auto new_size = ts2.cut(t1, t2, &new_epochs);
  ts2.set_epoch_vec(&new_epochs);
  assert((std::size_t)new_size == new_epochs.size() &&
         new_epochs.size() == ts2.data_points_vec().size());
  assert(new_epochs[0] >= t1 && new_epochs[new_size - 1] < t2);
  for (int i = 1; i < 12; i++) {
    auto dt1 = ngpt::datetime<ngpt::milliseconds>{
        ngpt::year(2019), ngpt::month(i), ngpt::day_of_month(1)};
    auto dt2 = ngpt::datetime<ngpt::milliseconds>{
        ngpt::year(2019), ngpt::month(i + 1), ngpt::day_of_month(1)};
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

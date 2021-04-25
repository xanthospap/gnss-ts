#include "eigen3/Eigen/Core"
#include "kalman.hpp"
#include "timeseries.hpp"
#include <ggdatetime/datetime_write.hpp>
#include <ggdatetime/dtcalendar.hpp>
#include <random>

std::default_random_engine e;
std::uniform_real_distribution<> dis(0, 1);

void dummy_print(const ngpt::timeseries &ts) {
  std::size_t idx = 0;
  while (idx < ts.size()) {
    std::cout << "\"" << ngpt::strftime_ymd_hmfs(ts.epochs_vec()[idx]) << "\" "
              << ts.data_points_vec()[idx].value() << "\n";
    ++idx;
  }
  return;
}

int main() {

  // vector of epochs
  std::vector<ngpt::datetime<ngpt::milliseconds>> epochs;
  ngpt::datetime<ngpt::milliseconds> t{ngpt::year(2016), ngpt::day_of_year(1)};
  ngpt::datetime_interval<ngpt::milliseconds> dt{ngpt::modified_julian_day(1),
                                                 ngpt::milliseconds(56000)};
  while (t < ngpt::datetime<ngpt::milliseconds>{ngpt::year(2020),
                                                ngpt::day_of_year(365)}) {
    epochs.emplace_back(t);
    t += dt;
  }

  // make a simple model
  ngpt::ts_model model;
  model.add_jump(ngpt::datetime<ngpt::milliseconds>{ngpt::year(2017),
                                                    ngpt::day_of_year(65)},
                 .123456e0);
  model.add_jump(ngpt::datetime<ngpt::milliseconds>{ngpt::year(2019),
                                                    ngpt::day_of_year(165)},
                 -0.423456e0);
  model.x0() = 1e-3;
  model.vx() = 0.0234567e0;

  // assign central epoch in model
  model.reference_epoch() = epochs[0];

  // create a time-series instance using the epoch vector and the model, adding
  // some white noise with std_dev = 0.01
  ngpt::timeseries ts(model, &epochs, .01);

  // make a model estimation
  ngpt::ts_model estimated(model);
  // estimated.zero_out_params();
  assert(estimated.num_parameters() == 4);
  std::vector<double> deviations = {1e-3, 1e-2, 1e-1, 1e0, 10e0, 100e0, 1e3};
  for (auto dev : deviations) {
    auto x = ngpt::kalman(ts, estimated, dev);
    std::cout << "\nRMS: " << std::sqrt(x.dot(model.state_vector()));
  }

  // print
  // dummy_print(ts);

  std::cout << "\n";
  return 0;
}

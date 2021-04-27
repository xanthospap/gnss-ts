#include "eigen3/Eigen/Dense"
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
  ngpt::timeseries ts(model, &epochs, .04);

  // make a model estimation (we want to start with a zero state vector)
  ngpt::ts_model estimated(model);
  std::vector<double> deviations = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 10e0, 100e0, 1e3};
  Eigen::VectorXd x = model.state_vector();
  for (auto dev : deviations) {
      auto state = ngpt::kalman(ts, estimated, dev);
      std::cout<<"\nSigma0 = "<<dev;
      for (int i=0; i<x.rows(); i++) std::cout<<"\n"<<std::abs(x(i)-state(i));
      std::cout<<"\n\tRMS for sigma0="<<dev<<" is: "<<std::sqrt(x.dot(state));
  }

  // print
  // dummy_print(ts);

  std::cout << "\n";
  return 0;
}

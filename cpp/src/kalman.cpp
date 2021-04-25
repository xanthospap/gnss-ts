#include "kalman.hpp"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/QR"
#include "ggdatetime/dtcalendar.hpp"
#ifdef DEBUG
#include <cassert>
#endif

Eigen::VectorXd ngpt::kalman(const ngpt::timeseries &ts,
                             const ngpt::ts_model &mdl, double sigma0) {
  auto num_parameters = mdl.num_parameters();
  Eigen::MatrixXd F = Eigen::RowVectorXd(num_parameters);
  Eigen::MatrixXd x = mdl.state_vector();
  Eigen::MatrixXd S = Eigen::MatrixXd(1, 1);
  Eigen::MatrixXd V = Eigen::MatrixXd(1, 1);
  V(0, 0) = 1e-2;
  Eigen::MatrixXd K = Eigen::MatrixXd(num_parameters, 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(num_parameters, num_parameters);
  Eigen::MatrixXd P = mdl.weight_matrix(sigma0);
  std::cout << "\nInitial P is: \n" << P;

  double t0 = mdl.reference_epoch().as_mjd(), yk_hat;
  ngpt::data_point dtp;
  std::size_t idx = 0, iparam = 0;
  for (auto t = ts.epochs_vec().cbegin(); t != ts.epochs_vec().cend(); ++t) {
    assert(idx < ts.data_points_vec().size());
    dtp = ts.data_points_vec()[idx];
    if (dtp.flag().is_clean()) {
      iparam = 0;
      F(iparam++) = 1e0;
      F(iparam++) = (t->as_mjd() - t0) / 365.25e0;
      for (const auto &jit : mdl.jumps_vec()) {
        F(iparam++) = (*t >= jit.start()) ? 1e0 : 0e0;
      }
      yk_hat = dtp.value() - (F * x)(0, 0); // forecast
      S = F * P * F.transpose() + V;        // residual covariance
      K = P * F.transpose() * S.inverse();  // kalman gain
      x = x + K * yk_hat;
      P = (I - K * F) * P;
      assert(iparam <= num_parameters);
    }
    ++idx;
  }
  std::cout << "\nestimated state:" << x;
  return x;
}

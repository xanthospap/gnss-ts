#include "kalman.hpp"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "ggdatetime/dtcalendar.hpp"
#ifdef DEBUG
#include <cassert>
#endif

Eigen::VectorXd ngpt::kalman(const ngpt::timeseries &ts,
                             const ngpt::ts_model &mdl, double sigma0) {
  auto num_parameters = mdl.num_parameters();
  Eigen::MatrixXd F = Eigen::RowVectorXd(num_parameters);
  Eigen::VectorXd x = /*Eigen::VectorXd::Zero(num_parameters);*/mdl.state_vector();
  Eigen::MatrixXd S = Eigen::MatrixXd(1, 1);
  Eigen::MatrixXd V = Eigen::MatrixXd(1, 1);
  V(0, 0) = 1e-2;
  Eigen::MatrixXd z = Eigen::MatrixXd(1, 1);
  Eigen::MatrixXd K = Eigen::MatrixXd(num_parameters, 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(num_parameters, num_parameters);
  Eigen::MatrixXd P = /*Eigen::MatrixXd::Identity(num_parameters, num_parameters);*/mdl.covariance_matrix(sigma0);

  double t0 = mdl.reference_epoch().as_mjd();
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
      z(0,0) = dtp.value();
      z = z - F * x; // forecast
      S = F * P * F.transpose() + V;        // residual covariance
      K = P * F.transpose() * S.inverse();  // kalman gain
      x = x + K * z;
      P = (I - K * F) * P;
      assert(iparam <= num_parameters);
    }
    ++idx;
  }
  return x;
}

/*
Eigen::VectorXd ngpt::kalman2(const ngpt::timeseries &ts,
                             const ngpt::ts_model &mdl, double sigma0) {
  auto num_parameters = mdl.num_parameters();
  typedef Eigen::Matrix<double, 1, 1> Matrix1d;
  Matrix1d z;
  Eigen::MatrixXd H = Eigen::RowVectorXd(num_parameters);
  Eigen::VectorXd x = mdl.state_vector();
  Eigen::VectorXd xp = x;
  Eigen::MatrixXd S = Eigen::MatrixXd(1, 1);
  Eigen::MatrixXd R = Eigen::MatrixXd(1, 1);
  Q(0, 0) = 1e-2;
  Eigen::MatrixXd K = Eigen::MatrixXd(num_parameters, 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(num_parameters, num_parameters);
  Eigen::MatrixXd P = mdl.weight_matrix(sigma0);
  // std::cout << "\nInitial P is: \n" << P;

  double t0 = mdl.reference_epoch().as_mjd();
  ngpt::data_point dtp;
  std::size_t idx = 0, iparam = 0;
  for (auto t = ts.epochs_vec().cbegin(); t != ts.epochs_vec().cend(); ++t) {
    assert(idx < ts.data_points_vec().size());
    dtp = ts.data_points_vec()[idx];
    z(0,0)= dtp.value();
    if (dtp.flag().is_clean()) {
      iparam = 0;
      H(iparam++) = 1e0;
      H(iparam++) = (t->as_mjd() - t0) / 365.25e0;
      for (const auto &jit : mdl.jumps_vec()) {
        H(iparam++) = (*t >= jit.start()) ? 1e0 : 0e0;
      }
      xp = H * x;
      P = H * P * H.transpose() + Q;
      K = P * H * (H*P*H.transpose()).inverse();
      x = xp + K*(z-H*xp);
      P = (I-K*H)*P;
      assert(iparam <= num_parameters);
    }
    ++idx;
  }
  // std::cout << "\nestimated state:" << x;
  return x;
}
*/

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/QR"
#include "ggdatetime/dtcalendar.hpp"

void kalman(const ngpt::timeseries& ts,c onst ngpt::model& mdl) {
  auto num_parameters = mdl.num_parameters();
  Eigen::MatrixXd F = Eigen::VectorXd(num_parameters);
  Eigen::MatrixXd x = Eigen::VectorXd(num_parameters);

  datetime<milliseconds> et; 
  double t0=mdl.reference_epoch().as_mjd();
  data_point dtp;
  std::size_t idx = 0, iparam=0
  for (auto it = m_epochs->begin(); it != m_epochs->end(); ++it) {
    et = *it;
    dtp = m_data[idx];
    iparam=0;
    if (dtp.flag().is_clean()) {
      F(iparam++) = 1e0;
      F(iparam++) = t0 - et.as_mjd();
      for (auto jit = m_jumps.begin(); jit != m_jumps.end(); ++m_jumps) {
        F(iparam++) = (et>=jit->start()) ? 1e0 : 0e0;
      }

    }
  }
  
}

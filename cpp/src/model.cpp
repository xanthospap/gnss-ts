#include "model.hpp"
#ifdef DEBUG
#include <ggdatetime/datetime_write.hpp>
#endif

std::size_t ngpt::ts_model::num_parameters() const noexcept {
  std::size_t parameters = 2;
  parameters += m_jumps.size();
  parameters += m_harmonics.size() * 2;
  parameters += m_vel_changes.size();
  return parameters;
}

std::size_t ngpt::ts_model::fill_data_vec(
    const std::vector<ngpt::datetime<ngpt::milliseconds>> *epochs,
    std::vector<ngpt::data_point> &data_vec) const noexcept {
  double t0 = m_reference_epoch.as_mjd();
  for (const auto &eph : *epochs) {
    double v = m_x0, t = (eph.as_mjd() - t0) / 365.25e0;
    v += m_vx * t;
    for (const auto &j : m_jumps) {
      v += (eph >= j.start()) ? j.value() : 0e0;
    }
    data_vec.emplace_back(v, 0e0);
  }
  return data_vec.size();
}

void ngpt::ts_model::zero_out_params() noexcept {
  m_x0 = m_vx = 0e0;
  for (auto &jump : m_jumps)
    jump.value() = 0e0;
  for (auto &harm : m_harmonics) {
    harm.in_phase() = harm.out_of_phase() = 0e0;
  }
  for (auto &velc : m_vel_changes)
    velc.value() = 0e0;
  return;
}

Eigen::MatrixXd ngpt::ts_model::weight_matrix(double sigma0) const noexcept {
  auto num_parameters = this->num_parameters();
  Eigen::MatrixXd P = Eigen::MatrixXd::Identity(num_parameters, num_parameters);
  std::size_t idx = 0;
  P(idx, idx) = sigma0 / m_x0_stddev;
  ++idx;
  P(idx, idx) = sigma0 / m_vx_stddev;
  ++idx;
  for (auto &jit : m_jumps) {
    P(idx, idx) = sigma0 / jit.stddev();
    ++idx;
  }
  return P;
}

Eigen::VectorXd ngpt::ts_model::state_vector() const noexcept {
  auto num_parameters = this->num_parameters();
  auto x0 = Eigen::VectorXd(num_parameters);
  std::size_t idx = 0;
  x0(idx++) = m_x0;
  x0(idx++) = m_vx;
  for (auto &jit : m_jumps) {
    x0(idx++) = jit.value();
  }
  return x0;
}
